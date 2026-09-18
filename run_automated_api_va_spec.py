"""Run the automated criteria through the API path and write their VA-Spec evidence lines.

The API path, not the CLI's. It parses the VCF with parse_request_vcf() and evaluates with
run_selected_criteria() - the two functions api/main.py calls for a
/v1/get_evidence_line_by_target_criteria request - in process rather than over HTTP, so no
server has to be running and the code that answers is the same either way.

Only the 16 automated criteria are asked for. The other 12 need PubMed and the LLM through
MCP, and evaluate_selected_criteria() accepts mcp=None precisely for a request like this one -
asking for a literature code without them raises rather than quietly returning UNKNOWN.

Policy and providers come from config/demo-rules.json, which is where the API reads them and
where they can be changed. Needs the network on a cold cache.
"""
import argparse
import asyncio
import json
from dataclasses import replace
from pathlib import Path

from acmg_pipeline.automated_engine import CRITERIA as AUTOMATED_CRITERIA
from acmg_pipeline.pipeline_interface import (
    load_automated_config, parse_request_vcf, run_selected_criteria,
)
from acmg_pipeline.providers.ensembl import EnsemblIdentityProvider
from acmg_pipeline.providers.http import CachedHttpClient

ROOT = Path(__file__).resolve().parent
RULES_PATH = ROOT / "config" / "demo-rules.json"
DEFAULT_INPUT = ROOT / "democase" / "case1_api_input_case1-var1.json"


def check_config():
    """The providers are the server's own settings; this only reports what they are.

    load_automated_config() reads config/demo-rules.json and says in its own docstring that
    callers cannot override it, so enabling a provider is a change to that file rather than
    an argument here - which is the point: a run through this script sees exactly what a run
    through the API sees.
    """
    config = load_automated_config()
    enabled = sorted(key for key, value in config.items()
                     if key.startswith("with_") and value)
    print(f"[api] providers enabled in config/demo-rules.json: {', '.join(enabled)}")
    return config


def requests_from(path):
    """The API request bodies to evaluate: one per data line of the request's own VCF.

    democase/*_api_input_*.json is an actual request body - {"vcf", "clinical_note"} - so this
    reads what a caller would post. The VCF carries several variants; the API answers about
    one, so each data line is sent on its own with the header it needs.
    """
    request = json.loads(path.read_text(encoding="utf-8"))
    header, bodies = [], []
    for line in request["vcf"].splitlines():
        if line.startswith("#"):
            header.append(line)
        elif line.strip():
            bodies.append(line)
    return ["\n".join([*header, body]) + "\n" for body in bodies]


def resolved(variant, provider):
    """The same variant with REF/ALT filled in, when the request left ALT as ".".

    The demo request bodies use the project's synthetic convention - ALT "." with the real
    identity in TRANSCRIPT/HGVSC - and the API path needs concrete alleles before it does
    anything at all. Resolving that is what the CLI's prepare step does and what the API has
    no step for, so it is done here, against the same Ensembl provider and release the
    evidence comes from. Everything after this is the API's own code.
    """
    if variant.alt != ".":
        return variant
    _candidate, annotation, _predictions = provider.map_record_with_evidence(
        {"identity": dict(variant.info)})
    assembly, chrom, pos, ref, alt = annotation["variant_key"].split(":")
    return replace(variant, chrom=chrom, pos=int(pos), ref=ref, alt=alt)


async def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "va_spec_output" / "automated_api_run.json")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    config = check_config()
    criteria = tuple(sorted(AUTOMATED_CRITERIA))
    bodies = requests_from(args.input)[:args.limit]
    print(f"[api] {len(bodies)} requests x {len(criteria)} automated criteria")
    client = CachedHttpClient(config.get("evidence_cache_dir")
                              or str(ROOT / "cache" / "erepo_automated_evidence"))
    provider = EnsemblIdentityProvider(
        client, config.get("ensembl_release")
        or EnsemblIdentityProvider.current_release(client))

    records = []
    for body in bodies:
        variant = resolved(parse_request_vcf(body), provider)
        # mcp/erepo_client stay None: an automated-only request never reaches PubMed or the
        # LLM, and run_selected_criteria() raises rather than quietly returning UNKNOWN if a
        # literature code is asked for without them.
        lines = await run_selected_criteria(
            variant, "", criteria, mcp=None, erepo_client=None,
        )
        records.append({
            "record_id": variant.id,
            "variant": {"assembly": "GRCh38", "chrom": variant.chrom, "pos": variant.pos,
                        "ref": variant.ref, "alt": variant.alt},
            "evidence_lines": lines,
        })
        applied = [code for code, line in lines.items()
                   if line.get("directionOfEvidenceProvided") in ("supports", "disputes")]
        print(f"  {variant.id} {variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt}"
              f" -> {len(lines)} lines, applied: {', '.join(applied) or 'none'}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps({"schema_version": "1.0", "source": "api path, automated criteria only",
                    "criteria": criteria, "records": records},
                   ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    asyncio.run(main())
