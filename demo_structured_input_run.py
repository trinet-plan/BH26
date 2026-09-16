"""
demo_structured_input_run.py

Demonstrates acmg_pipeline.pipeline.judge_variant_from_structured_input():
builds an ApiCaseInput the way a real API caller would - a plain dict with
"vcf"/"clinical_note" keys, exactly what json.loads(request_body) would
produce, never a file (ApiCaseInput.from_json_file() is a demo-fixture
convenience only, not the production path - see api_input.py) - then runs
it through the literature judgment pipeline.

Uses the real democase Case 3 MYH7 c.2155C>T variant (democase/case3_
variants_v2.vcf's actual case3-var1 row + INFO declarations) - this
variant is ClinGen 3-star reviewed with real ERepo PMIDs already in this
project's ground truth work, and its full text is already in cache/
pubmed_fulltext/ from the 227-pair validation run, so this demo runs fast
(cache hits) rather than needing fresh PubMed MCP fetches.

Deliberately literature-path only (PS3/BS3/PS4/PP1/BS4) - see pipeline.
judge_variant_from_structured_input()'s docstring for why clinical_note
isn't touched by this call at all.
"""

import asyncio
import sys
from contextlib import AsyncExitStack
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.api_input import ApiCaseInput
from acmg_pipeline.fulltext_cache import DiskBackedFullTextCache
from acmg_pipeline.gate import ERepoClient
import acmg_pipeline.pipeline as pl

# The embedded single-variant VCF text, exactly as democase/case3_variants_v2.vcf
# declares it (header meta-lines trimmed to what parse_vcf() actually reads) -
# this is what a real caller's JSON body's "vcf" field would contain.
_CASE3_VAR1_VCF = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description="RefSeq transcript">
##INFO=<ID=HGVSC,Number=1,Type=String,Description="HGVS coding">
##INFO=<ID=HGVSP,Number=1,Type=String,Description="HGVS protein">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical significance">
##INFO=<ID=CLNVARIATIONID,Number=1,Type=String,Description="ClinVar accession">
##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="Overall gnomAD AF">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
14\t23425971\tcase3-var1\tG\tA\t99\tPASS\tGENE=MYH7;TRANSCRIPT=NM_000257.4;HGVSC=c.2155C>T;HGVSP=p.(Arg719Trp);CLNSIG=Pathogenic;CLNVARIATIONID=VCV000014104;GNOMAD_AF=0.000005
"""

# The api_input.json envelope this represents - literally what an API
# caller's request body would deserialize into. clinical_note is included
# for completeness (ApiCaseInput carries it) but unused by this literature-
# only demo.
_API_REQUEST_BODY = {
    "vcf": _CASE3_VAR1_VCF,
    "clinical_note": Path("democase/case3_clinical_note_v1.txt").read_text(encoding="utf-8")
    if Path("democase/case3_clinical_note_v1.txt").exists() else "",
}


async def main() -> None:
    case_input = ApiCaseInput.from_dict(_API_REQUEST_BODY)  # dict in, no file involved
    pl.show(f"[demo] ApiCaseInput built from a dict (not a file): "
             f"vcf={len(case_input.vcf)} chars, clinical_note={len(case_input.clinical_note)} chars")

    full_text_cache = DiskBackedFullTextCache("cache/pubmed_fulltext")

    async with AsyncExitStack() as stack:
        mcp = await pl.connect_pubmed(stack)
        pl.show("[MCP] Connected to PubMed")
        erepo_client = ERepoClient()

        results = await pl.judge_variant_from_structured_input(
            case_input, mcp, erepo_client, full_text_cache=full_text_cache,
        )

        pl.show(f"\n{'='*70}\n[Results]\n{'='*70}")
        for criterion, aggregated in results.items():
            pl.show(f"  {criterion}: {aggregated.aggregated_direction.value}")


if __name__ == "__main__":
    asyncio.run(main())
