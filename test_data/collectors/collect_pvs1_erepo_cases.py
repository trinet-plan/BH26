"""Rebuild test_data/pvs1_erepo_cases.json from the live ClinGen ERepo API.

Queries each gene's classifications, keeps every interpretation carrying a PVS1
evidence code, and records what the API actually returns. Nothing here is
inferred except `variant_type_from_hgvs`, which is a mechanical read of the HGVS
string and is labelled as such.
"""
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
API = "https://erepo.clinicalgenome.org/evrepo/api/classifications"
GENES = ["PAH", "MYBPC3", "MYH7", "BRCA1", "BRCA2", "ATM", "PTEN", "CDH1", "RUNX1",
         "TP53", "LDLR", "GAA", "USH2A", "RPE65", "MYOC", "PALB2", "VHL", "DICER1",
         "KCNQ1", "HNF1A", "FBN1"]


def variant_type(hgvs_c, hgvs_p):
    """A mechanical read of the HGVS string - not a VEP consequence."""
    if re.search(r"c\.-?\d+[+-][12](?![0-9])", hgvs_c):
        return "canonical_splice"
    if re.fullmatch(r".*c\.[123][ACGT]>[ACGT]", hgvs_c) or re.search(r"p\.Met1(?![0-9])", hgvs_p):
        return "start_lost"
    if "fs" in hgvs_p:
        return "frameshift"
    if re.search(r"Ter|\*", hgvs_p):
        return "stop_gained"
    if re.search(r"del|dup|ins", hgvs_c):
        return "indel_unclassified"
    return "other"


def pick(hgvs_list, needle):
    return next((h for h in hgvs_list if needle in h), "")


def main():
    retrieved_at = datetime.now(timezone.utc).isoformat()
    entries = []
    for gene in GENES:
        response = requests.get(API, params={"gene": gene}, timeout=30)
        response.raise_for_status()
        for interp in response.json().get("variantInterpretations", []):
            codes = [
                (code.get("label"), code.get("status"))
                for guideline in interp.get("guidelines", [])
                for agent in guideline.get("agents", [])
                for code in agent.get("evidenceCodes", [])
            ]
            pvs1 = [(label, status) for label, status in codes
                    if label and label.upper().startswith("PVS1")]
            if not pvs1:
                continue
            hgvs = interp.get("hgvs", [])
            titled = pick(hgvs, "(")
            hgvs_c = pick(hgvs, ":c.")
            hgvs_p = titled if "p." in titled else ""
            label, status = pvs1[0]
            entries.append({
                "gene": gene,
                "hgvs_c": hgvs_c,
                "preferred_variant_title": titled,
                "caid": interp.get("caid"),
                "erepo_uuid": interp.get("uuid"),
                "erepo_url": interp.get("@id"),
                "condition": (interp.get("condition") or {}).get("label"),
                "pvs1_label": label,
                "pvs1_status": status,
                "classification": (interp.get("guidelines") or [{}])[0].get("outcome", {}).get("label"),
                "variant_type_from_hgvs": variant_type(hgvs_c, hgvs_p),
            })
        print(f"{gene}: {sum(1 for e in entries if e['gene'] == gene)} with a PVS1 code", flush=True)

    entries.sort(key=lambda e: (e["gene"], e["hgvs_c"]))
    document = {
        "schema_version": "1.0",
        "registry_status": "DRAFT",
        "purpose": (
            "Real ClinGen ERepo PVS1 outcomes, collected so that decision-tree coverage can be "
            "built from actual expert-panel calls instead of invented ones. These are expected "
            "OUTCOMES for comparison only. They MUST NOT be supplied to the evaluator as input: "
            "an expert-panel label that becomes its own supporting evidence is the circularity "
            "tests/test_clingen_positive.py checks for."
        ),
        "limitations": [
            "The ERepo API exposes a PVS1 evidence code as {label, status} only. It does not "
            "expose the decision-tree node values PVS1 needs - NMD prediction, exon relevance, "
            "critical-region disruption, protein loss fraction, splice or initiation assessment - "
            "so no fixture can be derived from this file alone.",
            "variant_type_from_hgvs is read off the HGVS string, not from a VEP consequence, and "
            "is a starting point for selection rather than an annotation.",
            "Genes were chosen by hand; this is not an exhaustive ERepo export.",
        ],
        "source": "ClinGen Evidence Repository API",
        "source_url": API,
        "retrieved_at": retrieved_at,
        "genes_queried": GENES,
        "entries": entries,
    }
    out = ROOT / "test_data" / "fetched_data" / "pvs1_erepo_cases.json"
    out.write_text(json.dumps(document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"\n{len(entries)} entries -> {out}")


main()
