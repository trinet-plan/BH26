# GA4GH VA-Spec feedback: gaps found building an ACMG/AMP pipeline

Compiled 2026-09-18, while removing duplicate/custom fields from this
project's VA-Spec `EvidenceLine` output (see `doc/curator_hints_unification_proposal_ja.md`
for the internal side of that cleanup). Each item below was found by trying
to express something this project's ACMG/AMP 2015 evaluation pipeline
genuinely needs using only standard VA-Spec fields, and hitting a real
limit - confirmed directly against the installed Python reference
implementation, not just read from the schema docs.

**Package**: `ga4gh.va_spec` (Python reference implementation), v0.5.0a4
(repo: `ga4gh/va-spec-python`)

Issues 1-3 are field *constraints* (a standard field exists but can't hold
what we need); issues 4-5 are field *additions* (no standard field exists
for the concept at all).

---

## Issue 1: `Method.reportedIn` cannot represent multiple contributing rule documents

**Affected model**: `ga4gh.va_spec.base.core.Method.reportedIn`

**Problem**: `Method.reportedIn` is typed as `Union[Document, iriReference, None]` — it accepts exactly one document. In practice, a single evaluation method is often built on more than one published source. For example, our ACMG/AMP PVS1 assessment method cites three documents at once:
- Richards et al. 2015 (the base ACMG/AMP guideline)
- the ClinGen SVI PVS1 recommendation (2018)
- the ClinGen SVI splicing update (2023)

There's no way to attach all three via the standard `specifiedBy.reportedIn` field. We worked around this with a custom extension (`rulesUsed`, an array of document citations) alongside `reportedIn`, which duplicates whichever one document we did put in `reportedIn`.

**Suggested fix**: Allow `Method.reportedIn` to accept `Document | list[Document] | iriReference | None`, mirroring how other "one-or-many" GA4GH fields are typically modeled.

**Reproduction**:
```python
from ga4gh.va_spec.base.core import Method
Method(methodType="PVS1", name="test", reportedIn=[doc1, doc2, doc3])  # fails - only one Document accepted
```

---

## Issue 2: `EvidenceLine.hasEvidenceItems` has no generic/open evidence-item type

**Affected model**: `ga4gh.va_spec.base.core.EvidenceLine.hasEvidenceItems`

**Problem**: `hasEvidenceItems` is typed as `list[Union[CohortAlleleFrequencyStudyResult, ExperimentalVariantFunctionalImpactStudyResult, Statement, EvidenceLine, iriReference]]`. The two `StudyResult` subtypes are useful, standardized profiles, but they only cover population-frequency and functional-assay evidence. Many real evaluation pipelines produce OTHER kinds of structured supporting facts that aren't judgments or statements themselves — e.g., a UniProt protein-domain lookup, or a ClinVar text-mining search result — and none of the five allowed member types fit:
- The two `StudyResult` subtypes reject an unrecognized `type` value (confirmed: passing `{"type": "StudyResult", ...}` raises a discriminator-tag validation error — there is no generic/open `StudyResult` base type actually usable here).
- A bare `iriReference` is schema-valid but practically useless without a separate, out-of-band catalog to resolve it against (which most consumers of a single `EvidenceLine` document won't have).
- The only way we found to legitimately attach open-ended structured facts is to fabricate a synthetic nested `EvidenceLine` and stash the facts in ITS OWN `extensions` — a workaround, not a first-class fit, since a "fact" is not itself an evidence line with a direction.

**Suggested fix**: Either (a) add a generic/open `StudyResult` variant to the union (freeform `type`, `id`, `name`, and an `extensions`-style bag for domain-specific facts), or (b) document the nested-`EvidenceLine` pattern as the officially sanctioned way to attach non-standardized structured evidence, so implementers aren't left guessing.

**Reproduction**:
```python
from ga4gh.va_spec.base.core import EvidenceLine, Method
EvidenceLine(
    id="x", directionOfEvidenceProvided="neutral",
    hasEvidenceItems=[{"id": "x", "type": "StudyResult", "name": "n"}],
    specifiedBy=Method(methodType="PM1", name="test"),
)  # raises: Input tag 'StudyResult' found using 'type' does not match any of the expected tags
```

---

## Issue 3: `directionOfEvidenceProvided` cannot distinguish "not evaluated" from a genuine neutral finding

**Affected model**: `ga4gh.va_spec.base.core.Direction` / `EvidenceLine.directionOfEvidenceProvided`

**Problem**: `Direction` has exactly three values: `supports`, `neutral`, `disputes`. This works well when an evaluation genuinely concluded "neither" — but many real evaluation pipelines (ACMG/AMP criteria being our case) have a THIRD distinct outcome that isn't "neutral evidence," it's "no evaluation was possible at all" (missing input data, unimplemented logic, inconclusive literature). Both cases currently have nowhere to go but `neutral`, so a consumer reading only the standard field cannot tell:

- "This criterion was evaluated and the evidence is genuinely neutral/inapplicable" (e.g., a criterion that requires functional-assay data, and the assay came back ambiguous)
- "This criterion was never meaningfully evaluated" (e.g., missing transcript annotation, an unimplemented rule, or a literature search that found nothing usable)

We worked around this by adding a custom extension field (`bh26AssessmentDetails.status`, one of `met`/`not_met`/`unknown`) carrying the real state, with a separate `direction` field (our own, not to be confused with the standard one) preserving the pre-collapse value (e.g. `"none"`) for the `not_met`/`unknown` case. This is exactly the kind of information a generic Statement/EvidenceLine schema should be able to carry without a project-specific extension, since "was this actually evaluated" is a universal question for any evidence-line consumer, not something specific to ACMG/AMP.

**Suggested fix**: Either (a) add a fourth `Direction` value (e.g. `not_evaluated`), or (b) add a separate, optional boolean/enum field to `EvidenceLine` (e.g. `evaluated: bool` or `assessmentOutcome: enum{met, not_met, not_evaluated}`) that a consumer can check before trusting `directionOfEvidenceProvided`'s three-way read.

**Reproduction / illustration**:
```python
from ga4gh.va_spec.base.core import EvidenceLine, Method

# Case A: assay ran, genuinely inconclusive
EvidenceLine(id="a", directionOfEvidenceProvided="neutral", specifiedBy=Method(methodType="PS3", name="x"))

# Case B: never evaluated at all (e.g. missing transcript data)
EvidenceLine(id="b", directionOfEvidenceProvided="neutral", specifiedBy=Method(methodType="PS3", name="x"))

# Both serialize identically on the standard field - a consumer has no way
# to distinguish A from B without a project-specific extension.
```

---

## Issue 4: No standard field for a step-by-step decision/reasoning trace on an EvidenceLine

**Affected model**: `ga4gh.va_spec.base.core.EvidenceLine`

**Problem**: `EvidenceLine`'s standard fields (`description`, `specifiedBy`, `directionOfEvidenceProvided`, `strengthOfEvidenceProvided`, `evidenceOutcome`, `hasEvidenceItems`) capture the *inputs* and the *conclusion* of an evaluation, but there is no field for the *reasoning path* an automated/rule-based `Method` walked through to get from one to the other. This is a generic need for any rule-based or decision-tree-driven evaluation, not specific to ACMG/AMP: our PVS1 evaluator, for example, produces an ordered list of decision-tree nodes (node id, the rule it applied, the intermediate result, which evidence items it consulted) that a reviewer or a debugging tool needs in order to audit *why* the line reached its stated outcome, beyond what `description`'s free text can reliably carry in a structured, machine-readable way.

We currently carry this as a custom extension (`bh26AssessmentDetails.decisionTrace`, a list of `{node_id, name, result, value, evidence_ids, rule_id, rule_source, next_node, note}` objects). Since `description` is free text and `hasEvidenceItems` only lists what was consulted (not the order or logic that connected them to the outcome), there's no standard place for this today.

**Suggested fix**: Add an optional field to `EvidenceLine` (e.g. `reasoningSteps` or `decisionTrace`) holding an ordered list of generic step objects (id, description, outcome, and a reference to the evidence item(s) that step consulted) - loosely enough specified that any rule-engine-based Method can populate it without a domain-specific schema.

---

## Issue 5: No standard field for reviewer-facing caveats/warnings on an EvidenceLine

**Affected model**: `ga4gh.va_spec.base.core.EvidenceLine` (or its `InformationEntity` ancestor)

**Problem**: An automated or LLM-assisted evaluation frequently needs to flag something a human reviewer should check before trusting the line's stated outcome - e.g. "this used a default/unreviewed policy threshold," "only one paper contributed a supporting judgment," "the disease mechanism has not been curator-confirmed." This is conceptually distinct from the outcome itself (`directionOfEvidenceProvided`/`evidenceOutcome`) and from the evidence used (`hasEvidenceItems`) - it's metadata about the *reliability/completeness* of the evaluation. There's no standard field for it today, and this need is generic to any automated evidence-generation pipeline, not specific to our use case.

We currently carry this as a custom extension (`curatorHints`, a list of `{severity, category, message}` objects - severity ∈ {info, caution, warning}).

**Suggested fix**: Add an optional field to `EvidenceLine` (e.g. `caveats` or `reviewFlags`) holding a list of `{severity, message}` (or similar) objects, so a generic VA-Spec consumer can surface "things to double-check" without needing a project-specific extension name to look for.
