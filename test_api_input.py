import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.api_input import ApiCaseInput
from acmg_pipeline.clinical_note import ClinicalNoteExtraction

passed = 0
failed = 0


def check(label, cond):
    global passed, failed
    if cond:
        passed += 1
        print(f"  OK   {label}")
    else:
        failed += 1
        print(f"  FAIL {label}")


DEMOCASE_DIR = Path(__file__).resolve().parent / "democase"
CASE1_FILES = sorted(DEMOCASE_DIR.glob("case1_api_input_case1-*.json"))

print(f"[1] Loading {len(CASE1_FILES)} case1 api_input file(s)")
check("found at least the 7 known case1 files", len(CASE1_FILES) >= 7)

expected_genes = {
    "case1-noise2": "MYBPC3", "case1-noise3": "MYBPC3", "case1-noise4": "TNNI3",
    "case1-noise5": "LMNA", "case1-var1": "MYBPC3", "case1-var2": "KCNJ5", "case1-var3": "MYH7",
}

for path in CASE1_FILES:
    case = ApiCaseInput.from_json_file(path)
    parsed = case.parse_vcf()
    record = parsed.record
    check(f"{path.name}: record id matches filename ({record.id})", record.id in path.stem)
    expected_gene = expected_genes.get(record.id)
    check(f"{path.name}: GENE={record.info.get('GENE')} matches expected {expected_gene}",
          expected_gene is None or record.info.get("GENE") == expected_gene)
    check(f"{path.name}: meta has reference/source", "reference" in parsed.meta and "source" in parsed.meta)
    check(f"{path.name}: at least one ##INFO declaration parsed", len(parsed.info_defs) > 0)

print("\n[1b] parse_vcf() rejects a VCF that doesn't have exactly 1 data line")
no_records_vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
try:
    ApiCaseInput(vcf=no_records_vcf, clinical_note="").parse_vcf()
    check("parse_vcf() raised ValueError on 0 data lines", False)
except ValueError:
    check("parse_vcf() raised ValueError on 0 data lines", True)

print("\n[2] Type coercion via declared ##INFO defs")
noise2 = ApiCaseInput.from_json_file(DEMOCASE_DIR / "case1_api_input_case1-noise2.json")
r = noise2.parse_vcf().record
check("AM_PATHOGENICITY is a float (Type=Float, declared)", isinstance(r.info["AM_PATHOGENICITY"], float))
check("AG_IMPACT_PHRED is a float", isinstance(r.info["AG_IMPACT_PHRED"], float))
check("GENE is a str", isinstance(r.info["GENE"], str))
check("GNOMAD_AF (undeclared in ##INFO) still coerces to float via fallback",
      isinstance(r.info["GNOMAD_AF"], float) and r.info["GNOMAD_AF"] == 0.0130)

print("\n[3] Number='.' INFO fields split into lists")
var2 = ApiCaseInput.from_json_file(DEMOCASE_DIR / "case1_api_input_case1-var2.json")
r2 = var2.parse_vcf().record
check("ACMG_CODES (Number=.) is a list", isinstance(r2.info["ACMG_CODES"], list))
check("ACMG_CODES has 2 entries for var2", r2.info["ACMG_CODES"] == ["PM2_Moderate", "PP3_Supporting"])
check("TAVTIGIAN_POINTS (undeclared) coerces to int", r2.info["TAVTIGIAN_POINTS"] == 3)

print("\n[4] ALT='.' (var1, a deletion) does not crash and is preserved as-is")
var1 = ApiCaseInput.from_json_file(DEMOCASE_DIR / "case1_api_input_case1-var1.json")
r1 = var1.parse_vcf().record
check("var1 ALT is '.'", r1.alt == ".")
check("var1 HGVSC is c.278delA", r1.info["HGVSC"] == "c.278delA")

print("\n[5] ClinicalNoteExtraction.from_json (parsing an LLM's structured output)")
fake_llm_output = {
    "proband": {
        "phenotype": {
            "affected_status": True,
            "clinical_features": [{"label": "palpitations", "hpo_id": None}],
            "age": "63",
            "age_of_onset": None,
        },
        "genotype": {"variant_status": True, "zygosity": "heterozygous"},
    },
    "family": {
        "inheritance_pattern": None,
        "family_history_status": False,
        "relatives": [
            {
                "relationship": "son",
                "biological_relationship_confirmed": None,
                "affected_status": False,
                "clinical_features": [],
                "age": None,
                "age_of_onset": None,
                "variant_status": False,
                "zygosity": None,
            },
        ],
    },
    "de_novo": {
        "father_variant_status": None,
        "mother_variant_status": None,
        "father_phenotype": None,
        "mother_phenotype": None,
        "paternity_confirmed": None,
        "maternity_confirmed": None,
    },
}
extraction = ClinicalNoteExtraction.from_json(fake_llm_output)
check("proband.phenotype.affected_status round-trips", extraction.proband.phenotype.affected_status is True)
check("proband.phenotype.clinical_features round-trips with null hpo_id",
      extraction.proband.phenotype.clinical_features[0].label == "palpitations"
      and extraction.proband.phenotype.clinical_features[0].hpo_id is None)
check("proband.phenotype.age round-trips", extraction.proband.phenotype.age == "63")
check("proband.genotype.zygosity round-trips", extraction.proband.genotype.zygosity == "heterozygous")
check("family.family_history_status round-trips", extraction.family.family_history_status is False)
check("1 relative extracted", len(extraction.family.relatives) == 1)
check("relative (son) affected_status=False round-trips", extraction.family.relatives[0].affected_status is False)
check("relative (son) variant_status=False round-trips", extraction.family.relatives[0].variant_status is False)
check("de_novo defaults are all None when not stated", extraction.de_novo.paternity_confirmed is None)

print(f"\n{'='*70}\n{passed} passed, {failed} failed\n{'='*70}")
sys.exit(1 if failed else 0)
