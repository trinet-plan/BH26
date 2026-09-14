# ACMG Criteria Automation: Tool Comparison

**version: v1**

Prepared: 2026-09-14
Purpose: Reference material for BH26 ExpertBoard. Compares how six existing tools (InterVar, AutoPVS1, AutoGVP, BIAS-2015, AcmGENTIC, Diablo_ACMG) actually derive each ACMG/AMP criterion, and confirms that none of them output in the GA4GH VA-Spec schema natively.

---

## Criterion-by-Criterion Comparison

| Criterion | InterVar | AutoPVS1 | AutoGVP | BIAS-2015 | AcmGENTIC | Diablo_ACMG |
|---|---|---|---|---|---|---|
| **PVS1** (loss of function) | Coarse gene-level rules (only checks LOF intolerance; a known weakness) | **Dedicated decision tree** (sophisticated logic accounting for predicted NMD, exon position, critical domain overlap, etc.) | **Overrides InterVar's result with AutoPVS1's result** (adopts the more refined one) | Custom implementation via Nirvana/VEP annotations | N/A (PS3/BS3-only tool) | open-cravat based (details unclear) |
| **PM2/BA1/BS1/BS2** (frequency) | **Threshold-based** lookup against gnomAD/ExAC etc. | N/A | Inherits InterVar's result | References gnomAD etc. via Nirvana/VEP, custom thresholds | N/A | Similar frequency DB lookup |
| **PP3/BP4** (in-silico) | Consensus across multiple tools (threshold-based) | N/A | Inherits InterVar | **Independently integrates AlphaMissense etc.** (as supplementary Nirvana data) | N/A | Similar |
| **PM1** (functional domain) | **Simple lookup** against known UniProt domain annotations | N/A | Inherits InterVar | Custom implementation | N/A | Similar |
| **PM5/PS1** (known variant at same residue) | **Lookup** against ClinVar | N/A | Inherits InterVar (+ also separately reflects ClinVar's own values) | Uses **dedicated lookup files** (e.g. `PS4_clinvar_submitter_counts.tsv`) | N/A | Similar |
| **PS4** (case-control) | **Lookup against a pre-built database, GWASdb v2, only** (cannot capture new studies) | N/A | Inherits InterVar | Uses a **dedicated file tallying ClinVar submitter counts** (broader coverage than InterVar) | N/A | Similar |
| **PP5/BP6** (third-party assertions) | **Directly reflects ClinVar's own assertion** (no independent evaluation) | N/A | **Directly uses real ClinVar data** (star rating, submitter breakdown), with a conflict-resolution script | Similarly relies on ClinVar | N/A | Similar |
| **PS2, PS3, PM3, PM6, PP1, PP4** (require literature reading / pedigree data) | **Cannot be automated at all** (9-10 criteria left blank) | N/A | Likewise left blank (assumes manual input) | Likewise left blank; manual completion required via `user_classifiers` | **Only PS3/BS3 automated via LLM** (direction: high accuracy; strength: low accuracy) | Likewise left blank |

---

## Observations

1. **Even the same-named "PVS1" is judged with very different sophistication between InterVar and AutoPVS1** — which is exactly why AutoGVP is designed to override InterVar's result with AutoPVS1's when they disagree.
2. **PS4 shares the same underlying idea ("lookup against a pre-built database") across tools, but the actual source differs** (InterVar uses GWASdb; BIAS-2015 uses a tally of ClinVar submitter counts) — meaning the same variant can hit or miss depending on which tool is used.
3. **PP5/BP6 is handled by every tool the same way: simply trusting whatever ClinVar already asserts**, with no independent verification of its own.
4. **The six criteria requiring literature reading (PS2, PS3, PM3, PM6, PP1, PP4) are out of reach for every tool except AcmGENTIC** — which matches exactly the set of criteria we carved out as "Layer 3" in our own design, reinforcing that our choice of where to apply LLM support was well-targeted, as evidenced by the limitations of existing tools themselves.

---

## Related finding: none of these tools output GA4GH VA-Spec format

All six tools produce their own custom output format (TSV in most cases; AcmGENTIC produces a custom JSON with `decision`/`strength`/`narrative` fields). None natively emit GA4GH VA-Spec JSON. See `BH26_ExpertBoard_reference_material` for the fuller VA-Spec investigation, including a working code example that maps our own Evidence Record fields onto the VA-Spec `acmg_2015` Pydantic schema (`ga4gh.va_spec`), validated against real BH26 Demo Case data.
