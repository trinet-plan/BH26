"""
TogoMCP + PubMed MCP サンプル(複数MCPサーバー対応版)
【変更点】
  ・TogoMCPに加えて、PubMed MCPにも同時接続
  ・ツール名は "サーバー名__元のツール名" の形でOpenAIに渡し、
    衝突を避けつつ呼び出し時に正しいサーバーへルーティングする
【出力】
  画面    : 1周1行のサマリ (経過秒、トークン数、呼んだツール)
  ログ    : logs/run_YYYYmmdd_HHMMSS.log に全文 (クエリ・結果・最終回答)

【必要なもの】
  pip install mcp openai
"""

import asyncio
import json
import os
import time
from contextlib import AsyncExitStack
from datetime import datetime
from pathlib import Path

from mcp import ClientSession
from mcp.client.streamable_http import streamable_http_client
from openai import OpenAI

# ---------------------------------------------------------------------------
# 設定
# ---------------------------------------------------------------------------
ROOT_DIR = Path(__file__).resolve().parent


def _load_env_file(path: Path = ROOT_DIR / ".env") -> None:
    """.env があれば os.environ に読み込む(python-dotenv 不使用の簡易実装)。"""
    if not path.exists():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip())


_load_env_file()

VLLM_BASE_URL = os.environ.get("VLLM_BASE_URL", "")
VLLM_API_KEY = os.environ.get("VLLM_API_KEY", "")
if not VLLM_BASE_URL or not VLLM_API_KEY:
    raise RuntimeError(
        "VLLM_BASE_URL / VLLM_API_KEY が設定されていません。"
        "リポジトリ直下に .env を作成してください(.env.example を参照)。"
    )

# サーバーにロードされているモデルにあわせる
MODEL = "google/gemma-4-26B-A4B-it"
#MODEL = "RadixArk/Qwen3.8-27B-NVFP4"

# 接続するMCPサーバー一覧。ここに追加するだけで他のMCPサーバーも増やせる。
MCP_SERVERS = {
    "togomcp": "https://togomcp.rdfportal.org/mcp",
    "pubmed": "https://pubmed.mcp.claude.com/mcp",
}

MAX_TURNS = 15  # 暴走防止。これがないと同じツールを延々と呼び続けることがある

# これがないと「計画を書いて終わる」ことがある
SYSTEM_PROMPT = (
    "計画を文章で説明するのではなく、必要なツールを実際に呼び出して実行してください。"
    "最終的な答えが得られるまでツール呼び出しを続けてください。"
    "複数のツールソース(TogoMCP、PubMed)が使えるので、質問の内容に応じて"
    "適切な方を選んでください。両方必要な場合は両方使ってください。"
)

client = OpenAI(base_url=VLLM_BASE_URL, api_key=VLLM_API_KEY)


# ---------------------------------------------------------------------------
# ログ
# ---------------------------------------------------------------------------

LOG_DIR = Path("logs")
LOG_DIR.mkdir(exist_ok=True)
LOG_PATH = LOG_DIR / f"run_{datetime.now():%Y%m%d_%H%M%S}.log"


def log(text: str = ""):
    """ファイルにだけ書く"""
    with LOG_PATH.open("a", encoding="utf-8") as f:
        f.write(text + "\n")


def show(text: str):
    """画面とファイルの両方に書く"""
    print(text)
    log(text)


# ---------------------------------------------------------------------------
# MCPのツール定義 → OpenAIのfunction定義への変換
# ---------------------------------------------------------------------------
#
# 複数サーバーを束ねると、同じ名前のツールが別サーバーに存在する可能性がある
# (例: 両方に "search" というツールがある、等)。
# これを避けるため、OpenAIに渡す関数名は "サーバー名__元のツール名" にする。
# 呼び出し時はこのプレフィックスを見て、どのMCPセッションに投げるかを決める。

NAMESPACE_SEP = "__"


def namespaced_name(server_name: str, tool_name: str) -> str:
    return f"{server_name}{NAMESPACE_SEP}{tool_name}"


def split_namespaced_name(name: str) -> tuple[str, str]:
    server_name, _, tool_name = name.partition(NAMESPACE_SEP)
    return server_name, tool_name


def to_openai_tools(server_name: str, mcp_tools) -> list[dict]:
    return [
        {
            "type": "function",
            "function": {
                "name": namespaced_name(server_name, t.name),
                "description": f"[{server_name}] {t.description or ''}",
                "parameters": t.input_schema,
            },
        }
        for t in mcp_tools
    ]


def extract_text(result) -> str:
    """MCPのツール結果からテキストを取り出す"""
    return "\n".join(getattr(b, "text", str(b)) for b in result.content)


# ---------------------------------------------------------------------------
# 複数MCPサーバーへの同時接続
# ---------------------------------------------------------------------------


async def connect_all_mcp_servers(
    stack: AsyncExitStack,
) -> dict[str, ClientSession]:
    """MCP_SERVERSに列挙した全サーバーに接続し、名前 -> セッション の辞書を返す"""
    sessions: dict[str, ClientSession] = {}
    for server_name, url in MCP_SERVERS.items():
        read, write = await stack.enter_async_context(streamable_http_client(url))
        mcp = await stack.enter_async_context(ClientSession(read, write))
        await mcp.initialize()
        sessions[server_name] = mcp
    return sessions


# ---------------------------------------------------------------------------
# エージェントループ本体
# ---------------------------------------------------------------------------


async def ask(question: str) -> str:
    t_start = time.perf_counter()

    log(f"MODEL     : {MODEL}")
    log(f"QUESTION  : {question}")
    log(f"MCP SERVERS: {list(MCP_SERVERS.keys())}")
    log("=" * 70)

    async with AsyncExitStack() as stack:
        sessions = await connect_all_mcp_servers(stack)

        # 各サーバーからツール一覧を取得し、名前空間を付けて1つのリストにまとめる
        tools: list[dict] = []
        for server_name, mcp in sessions.items():
            server_tools = to_openai_tools(server_name, (await mcp.list_tools()).tools)
            tools.extend(server_tools)
            show(f"[MCP] {server_name}: ツール {len(server_tools)} 個を取得")
            log("  " + ", ".join(t["function"]["name"] for t in server_tools))

        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": question},
        ]

        for turn in range(1, MAX_TURNS + 1):
            t_turn = time.perf_counter()

            # --- (1) LLMに現在の履歴とツール一覧を渡して、次の一手を聞く ---
            res = client.chat.completions.create(
                model=MODEL,
                messages=messages,
                tools=tools,
            )
            msg = res.choices[0].message
            finish = res.choices[0].finish_reason

            dt = time.perf_counter() - t_turn
            elapsed = time.perf_counter() - t_start
            p, c = res.usage.prompt_tokens, res.usage.completion_tokens
            tps = c / dt if dt > 0 else 0

            head = (
                f"[{turn:2d}] {dt:6.1f}s (計{elapsed:6.1f}s) "
                f"p={p:<7} c={c:<6} {tps:5.1f}tok/s  {finish}"
            )

            log("\n" + "=" * 70)
            log(head)
            if getattr(msg, "reasoning", None):
                log(f"--- 思考 ---\n{msg.reasoning}")
            if msg.content:
                log(f"--- 本文 ---\n{msg.content}")

            # --- (2) LLMの発言を履歴に追加する ---
            #
            # これを忘れると、次のリクエストで「なぜtool結果が来ているのか」が
            # 分からなくなって壊れる。
            messages.append(
                {
                    "role": "assistant",
                    "content": msg.content,
                    "tool_calls": [
                        {
                            "id": c_.id,
                            "type": "function",
                            "function": {
                                "name": c_.function.name,
                                "arguments": c_.function.arguments,
                            },
                        }
                        for c_ in (msg.tool_calls or [])
                    ]
                    or None,
                }
            )

            # --- (3) ツール呼び出しでなければ、それが最終回答 ---
            if finish != "tool_calls":
                total = time.perf_counter() - t_start
                print(f"{head}  → 完了")
                print(f"\n{turn}周 / {total:.1f}秒 / ログ: {LOG_PATH}")
                log(f"\n完了: {turn}周 / {total:.1f}秒")
                return msg.content

            # --- (4) 要求されたツールを実際に呼ぶ(名前空間から対象サーバーを判定) ---
            names = []
            for call in msg.tool_calls:
                namespaced = call.function.name
                server_name, tool_name = split_namespaced_name(namespaced)
                args = json.loads(call.function.arguments)

                log(f"\n--- CALL {namespaced} (server={server_name}) ---")
                log(json.dumps(args, ensure_ascii=False, indent=2))

                t_tool = time.perf_counter()
                mcp = sessions.get(server_name)
                try:
                    if mcp is None:
                        raise RuntimeError(
                            f"未知のサーバー名 '{server_name}' "
                            f"(ツール名 '{namespaced}' の名前空間解析に失敗した可能性)"
                        )
                    content = extract_text(await mcp.call_tool(tool_name, args))
                except Exception as e:
                    # エラーもLLMに返す。そうすればLLM側で修正を試みられる
                    content = f"ERROR: {e}"
                tool_dt = time.perf_counter() - t_tool

                log(f"--- RESULT ({len(content)}字 / {tool_dt:.1f}s) ---")
                log(content)

                names.append(f"{namespaced}←{len(content)}字")

                # --- (5) 結果を role="tool" で履歴に追加 ---
                #
                # tool_call_id で「どの呼び出しへの答えか」を紐付ける。
                # 複数ツールを同時に呼ぶことがあるため必須。
                messages.append(
                    {
                        "role": "tool",
                        "tool_call_id": call.id,
                        "content": content,
                    }
                )

            print(f"{head}  → {', '.join(names)}")

            # ループの先頭に戻る。
            # messages が伸びた状態で、またLLMに聞き直す。

        total = time.perf_counter() - t_start
        print(f"\n{MAX_TURNS}周しても終わらなかった / {total:.1f}秒 / ログ: {LOG_PATH}")
        return f"{MAX_TURNS}周しても終わらなかった"


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    question = "UniProtでアルツハイマー病に関連するヒトタンパク質を探してください"

    answer = asyncio.run(ask(question))

    log("\n" + "=" * 70)
    log("--- 最終回答 ---")
    log(answer or "(なし)")

    print("\n" + "=" * 70)
    print(answer)
