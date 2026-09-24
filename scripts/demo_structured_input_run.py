"""
demo_structured_input_run.py

Demonstrates acmg_pipeline.pipeline.classify_variant_from_structured_input():
builds an ApiCaseInput the way a real API caller would - a plain dict with
"vcf"/"clinical_note" keys, exactly what json.loads(request_body) would
produce, never a file (ApiCaseInput.from_json_file() is a demo-fixture
convenience only, not the production path - see api_input.py) - runs the
implemented literature criteria (PS3/BS3/PS4 as of 2026-09-16 - PP1/BS4
were dropped from the default scope on that date, see acmg_pipeline.
classification.IMPLEMENTED_CODES), fills the other 25 with stubs.
stub_evidence() (their real logic is another team member's responsibility
- see acmg_pipeline/criteria/stubs.py), and classifies the result.

Two real democase scenarios, chosen to exercise both branches of judge_
variant_from_structured_input()'s ERepo lookup:

  1. MYH7 c.2155C>T (Case 3, case3-var1) - ClinGen 3-star reviewed, has real
     ERepo PMIDs (and its full text is already in cache/pubmed_fulltext/
     from the 227-pair validation run, so this branch runs fast).
  2. MYBPC3 c.278delA (Case 1, case1-var1) - a real but NOVEL variant (per
     democase/real_cases_groundtruth_integrated_v6_ja.md section 3: not
     registered in ClinVar/ClinGen at all). ERepoClient.lookup() returns
     zero PMIDs for it, so this exercises the "literature path finds
     nothing, everything falls to not_evaluated/stub" branch - the same
     gap discussed with the user on 2026-09-16 (no PubMed-search fallback
     exists yet for variants ERepo has never seen).

Each embedded VCF INFO string below deliberately OMITS the demo file's own
ACMG_CODES/TAVTIGIAN_POINTS fields - those are the demo data's own
precomputed answer key (see vcf_record.py's docstring), not real annotation
pipeline output, and must never be read as this project's own input.
"""

import asyncio
import sys
from contextlib import AsyncExitStack
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from acmg_pipeline.api_input import ApiCaseInput
from acmg_pipeline.fulltext_cache import DiskBackedFullTextCache
from acmg_pipeline.gate import ERepoClient
from acmg_pipeline.llm_cache import DiskBackedLLMCache
import acmg_pipeline.pipeline as pl

_VCF_HEADER = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description="RefSeq transcript">
##INFO=<ID=HGVSC,Number=1,Type=String,Description="HGVS coding">
##INFO=<ID=HGVSP,Number=1,Type=String,Description="HGVS protein">
##INFO=<ID=ZYGOSITY,Number=1,Type=String,Description="heterozygous/homozygous in proband">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical significance">
##INFO=<ID=CLNVARIATIONID,Number=1,Type=String,Description="ClinVar accession if registered">
##INFO=<ID=DISEASE_ASSOCIATION,Number=1,Type=String,Description="Primary disease association of the gene">
##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="Overall gnomAD AF">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

_SCENARIOS = {
    "MYH7 c.2155C>T (Case 3 - has ERepo PMIDs)": (
        _VCF_HEADER +
        "14\t23425971\tcase3-var1\tG\tA\t99\tPASS\t"
        "GENE=MYH7;TRANSCRIPT=NM_000257.4;HGVSC=c.2155C>T;HGVSP=p.(Arg719Trp);"
        "CLNSIG=Pathogenic;CLNVARIATIONID=VCV000014104;GNOMAD_AF=0.000005\n",
        "democase/case3_clinical_note_v1.txt",
    ),
    "MYBPC3 c.278delA (Case 1 - novel, zero ERepo PMIDs)": (
        _VCF_HEADER +
        "11\t47352561\tcase1-var1\tG\t.\t99\tPASS\t"
        "GENE=MYBPC3;TRANSCRIPT=NM_000256.3;HGVSC=c.278delA;HGVSP=p.(Lys93ArgfsTer3);"
        "ZYGOSITY=heterozygous;CLNSIG=Likely_Pathogenic;CLNVARIATIONID=not_registered_novel;"
        "DISEASE_ASSOCIATION=hypertrophic_cardiomyopathy\n",
        "democase/case1_clinical_note_v1.txt",
    ),
}


async def main() -> None:
    full_text_cache = DiskBackedFullTextCache("cache/pubmed_fulltext")
    llm_cache = DiskBackedLLMCache("cache/llm_judgments")

    async with AsyncExitStack() as stack:
        mcp = await pl.connect_pubmed(stack)
        pl.show("[MCP] Connected to PubMed")
        erepo_client = ERepoClient()

        for label, (vcf_text, note_path) in _SCENARIOS.items():
            pl.show(f"\n\n{'#'*70}\n# Scenario: {label}\n{'#'*70}")

            request_body = {
                "vcf": vcf_text,
                "clinical_note": Path(note_path).read_text(encoding="utf-8") if Path(note_path).exists() else "",
            }
            case_input = ApiCaseInput.from_dict(request_body)  # dict in, no file involved
            pl.show(f"[demo] ApiCaseInput built from a dict (not a file): "
                    f"vcf={len(case_input.vcf)} chars, clinical_note={len(case_input.clinical_note)} chars")

            await pl.classify_variant_from_structured_input(
                case_input, mcp, erepo_client, full_text_cache=full_text_cache, llm_cache=llm_cache,
            )


if __name__ == "__main__":
    asyncio.run(main())
