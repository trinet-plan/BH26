"""UniProt-derived first-pass answer for PM4/BP3's repeat-region question.

[What PM4/BP3 ask]
  An in-frame insertion, in-frame deletion, or stop-loss changes protein length. PM4
  supports pathogenicity unless the altered interval falls in a repetitive region
  without established function; BP3 supports a benign read exactly when it does. Both
  read the same `region` evidence category - see acmg_pipeline/criteria/regions.py -
  through `repetitive`/`functional_importance`/`nonfunctional_repeat`/
  `functional_review_complete`, fields that were previously never populated by any
  automated provider (the only other `region` producer, ClinVarHotspotProvider, only
  ever answers PM1's hotspot route - see providers/clinvar.py).

[The automated first pass]
  UniProt's own feature track answers this for many genes without new infrastructure,
  reusing the fetch/overlap helpers already built for PVS1's NF04/NF06 pass (see
  providers/uniprot_features.py and providers/protein_region.py, which this shares
  them with rather than reimplementing): a "Repeat" or "Compositional bias" feature
  over the altered interval (see REPETITIVE_FEATURE_TYPES below for why both count) is
  the repetitiveness this criterion means, and a Domain/Binding site/Active
  site/Motif/Coiled coil (or non-disordered Region) feature ALSO over that interval is
  the established-function counter-signal.

[Absence is evidence too, unlike protein_region.py's critical_region_disrupted]
  NF04/NF06 leave `critical_region_disrupted` unset when nothing overlaps, because "no
  named domain here" is not proof nothing important is there. Here, absence of a
  documented repetitive feature at the altered interval IS the signal repetitive=False
  records: PM4/BP3 ask about a documented repeat, not merely "any function unknown",
  so a resolved UniProt entry with no such feature there answers repetitive=False the
  same way protein_region.py's own region_biologically_relevant already treats an
  answered-but-empty UniProt entry as a real (not missing) result - see that module's
  docstring for the precedent. A record is still produced only once UniProt itself
  could be resolved; no accession or no feature list at all yields no record, same
  honest-gap convention as everywhere else in this project.

[Always flagged, never presented as a curator's own review]
  Every record here carries assessment_method="automated" with a named method/
  policy_version, exactly what require_boolean_fields's caller in regions.py accepts
  from a non-curator source - a curator's own reviewed region record still overrides
  this whenever one exists (curated_context selects on reviewed_or_automated(), not on
  who ran first).
"""

from __future__ import annotations

import hashlib

from acmg_pipeline.criteria.reference_links import gene_to_uniprot_accession
from acmg_pipeline.providers.uniprot_features import (
    CRITICAL_FEATURE_TYPES, fetch_features, feature_interval, is_disordered, overlaps,
)

METHOD = "uniprot_repeat_feature_overlap"
# UniProt's own vocabulary splits "repetitive" across two feature types: "Repeat" names a
# tandem structural repeat unit (ANK, WD40, LRR...), while a low-complexity/homopolymeric
# run (poly-Gln, poly-Pro, poly-Ala...) - the kind BP3's real-world examples are usually made
# of, e.g. FOXG1's proline-rich stretch and RPGR ORF15's compositionally biased tail - is
# instead typed "Compositional bias" (2026-09-21: found empirically after the 679-variant
# ground-truth run gave BP3 0% recall even on the FOXG1/RPGR cases it should have caught -
# both had a real UniProt feature over the altered interval, just not one this provider
# checked for). A Compositional bias run is, by UniProt's own definition, not a named
# structural domain, so unlike "Region" it is never also read as functional_importance.
REPETITIVE_FEATURE_TYPES = {"Repeat", "Compositional bias"}


class RegionRepeatProvider:
    """One `region` record answering PM4/BP3's repeat-region question, or none."""

    name = "UniProt repeat/domain overlap"

    def __init__(self, client, policy_version):
        if not policy_version:
            raise ValueError("policy_version must be recorded")
        self.client = client
        self.policy_version = policy_version

    def get_region(self, variant, gene, transcript, protein_id, start, end):
        if not (gene and transcript and isinstance(start, int) and start > 0):
            return []
        end = end if isinstance(end, int) and end >= start else start
        accession = gene_to_uniprot_accession(gene)
        if not accession:
            return []
        features, response = fetch_features(self.client, accession)
        if features is None:
            return []
        repetitive = False
        functional = False
        for feature in features:
            interval = feature_interval(feature)
            if interval is None or not overlaps(interval[0], interval[1], start, end):
                continue
            feature_type = feature.get("type")
            if feature_type in REPETITIVE_FEATURE_TYPES:
                repetitive = True
            elif feature_type in CRITICAL_FEATURE_TYPES and not (
                feature_type == "Region" and is_disordered(feature)
            ):
                functional = True
        digest = hashlib.sha256(
            f"{response['body_sha256']}:{start}-{end}".encode("utf-8")).hexdigest()
        return [{
            "category": "region",
            "variant_key": variant.key,
            "evidence_id": f"urn:sha256:{digest}:uniprot-repeat:{transcript}:{start}-{end}",
            "source": self.name,
            "source_version": accession,
            "retrieved_at": response["retrieved_at"],
            "quality_status": "PASS",
            "transcript": transcript,
            "protein_id": protein_id,
            "gene": gene,
            "start": start,
            "end": end,
            "region_type": "uniprot_repeat",
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": self.policy_version,
            "repetitive": repetitive,
            "functional_importance": functional,
            "nonfunctional_repeat": repetitive and not functional,
            "functional_review_complete": True,
            "response_sha256": response["body_sha256"],
            "policy_note": (
                "UniProt Repeat/Compositional-bias feature overlap over the altered protein "
                "interval, flagged as an automated first pass (assessment_method=automated), "
                "not a curator's own review. functional_importance is set when a Domain/"
                "Binding site/Active site/Motif/Coiled coil (or non-disordered Region) "
                "feature ALSO overlaps the same interval - see this provider's own module "
                "docstring for why no overlapping repetitive feature is read as "
                "repetitive=False rather than left unset."
            ),
        }]


__all__ = ["METHOD", "RegionRepeatProvider"]
