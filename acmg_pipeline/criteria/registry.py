"""
registry.py
Single source of truth for "does this project have a real judgment module
for ACMG code X, or only a stub" - ties acmg_pipeline.criteria.{ps3_bs3,
ps4, segregation} (real) together with acmg_pipeline.criteria.stubs
(everything else) via acmg_pipeline.classification's canonical code lists.

pipeline.py's ENGINE_BY_CRITERION (real JudgmentEngine instances - build_
prompt/from_json/finalize/aggregate, used to actually run the LLM) stays
there, not here: stub codes have no engine to run at all, so a registry
that only lists "which codes have engines" would be misleading about the
remaining seven. This module is about classification-time evidence lookup instead
- see get_criterion_evidence().
"""

from __future__ import annotations

from typing import Optional

from acmg_pipeline.classification import CriterionEvidence
from acmg_pipeline.constants import ALL_ACMG_CODES, IMPLEMENTED_CODES
from acmg_pipeline.criteria import stubs


def is_implemented(code: str) -> bool:
    if code not in ALL_ACMG_CODES:
        raise ValueError(f"{code!r} is not a recognized ACMG/AMP 2015 code")
    return code in IMPLEMENTED_CODES


def get_criterion_evidence(code: str, real_evidence: Optional[CriterionEvidence] = None) -> CriterionEvidence:
    """
    Returns the CriterionEvidence to actually use for `code` in a
    classify() call.

    - For an implemented literature or automated code, the caller must have
      already produced a real CriterionEvidence (e.g. via
      classification.from_aggregated_judgment() after running the LLM
      pipeline) and pass it as `real_evidence` - this function does not run
      anything itself.
    - For a stub code, `real_evidence` is ignored (there is nothing real to
      pass) and a NOT_EVALUATED placeholder is returned via
      acmg_pipeline.criteria.stubs.stub_evidence().

    This lets a caller building a full 28-code evidence set for classify()
    write one uniform loop over ALL_ACMG_CODES instead of branching on
    implemented-vs-stub itself:

        evidence = [
            get_criterion_evidence(code, my_results.get(code))
            for code in ALL_ACMG_CODES
        ]
    """
    if not is_implemented(code):
        return stubs.stub_evidence(code)
    if real_evidence is None:
        raise ValueError(
            f"{code!r} is an implemented code - it needs a real "
            "CriterionEvidence (e.g. from classification."
            "from_aggregated_judgment()), not None. Only stub codes get a "
            "placeholder automatically."
        )
    if real_evidence.code != code:
        raise ValueError(f"real_evidence.code={real_evidence.code!r} does not match requested code={code!r}")
    return real_evidence
