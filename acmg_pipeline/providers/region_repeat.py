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

[A second signal: sequence complexity, when features are silent]
  Not every real repeat is curated as a UniProt feature. RPGR's ORF15 region (the
  ground-truth source of most of this criterion's real BP3-met examples) has no
  Repeat/Compositional bias feature at all in UniProt's own data, despite its raw
  sequence there being a textbook low-complexity Glu/Gly-rich run. repetitive is
  therefore also set when the altered interval's own local sequence complexity - the
  Wootton-Federhen/SEG measure, see uniprot_features.sequence_complexity() - falls at
  or below SEG's own published trigger threshold (K1=2.2 bits). This never sets
  functional_importance; a low-complexity window is not evidence of an established
  domain either way.

[Isoform escalation, when a position falls outside the canonical sequence]
  A transcript's own protein numbering can run past the canonical UniProt entry's
  length when it corresponds to a longer alternatively-spliced isoform instead (RPGR-
  ORF15 is exactly this: residue ~1040 on NM_001034853.2, but the canonical entry
  Q92834 is only 1020 residues). uniprot_features.resolve_covering_entry() escalates
  to the first isoform whose own sequence is long enough to cover the altered
  interval, and the record's source_version/isoform_used name which accession the
  features and sequence actually came from.

[Absence is evidence too, unlike protein_region.py's critical_region_disrupted]
  NF04/NF06 leave `critical_region_disrupted` unset when nothing overlaps, because "no
  named domain here" is not proof nothing important is there. Here, absence of a
  documented repetitive feature at the altered interval IS the signal repetitive=False
  records: PM4/BP3 ask about a documented repeat, not merely "any function unknown",
  so a resolved UniProt entry with no such feature there answers repetitive=False the
  same way protein_region.py's own region_biologically_relevant already treats an
  answered-but-empty UniProt entry as a real (not missing) result - see that module's
  docstring for the precedent. A record is still produced only once an entry covering
  the altered interval could be resolved at all (no accession, or no entry long
  enough, yields no record); an entry that resolved but reports no `features` list
  still produces one, scored on sequence complexity alone, same honest-gap convention
  as everywhere else in this project.

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
    CRITICAL_FEATURE_TYPES, feature_interval, is_disordered, overlaps,
    resolve_covering_entry, sequence_complexity,
)

METHOD = "uniprot_repeat_feature_overlap"
# A second, sequence-derived repetitiveness signal alongside the feature-based one
# below - found necessary after RPGR's real UniProt entry (Q92834 canonical, and its
# ORF15 isoform Q92834-6) turned out to carry NO feature annotation at all over the
# ORF15 tail's own low-complexity Glu/Gly-rich stretch - the region is real (its raw
# sequence there, e.g. residues 1030-1045, reads "EGEGEEEEGEEEGEE"), UniProt's curators
# simply never tagged it as Repeat/Compositional bias. See uniprot_features.
# sequence_complexity()'s docstring for the Wootton-Federhen/SEG formula and citation.
#
# SEG_WINDOW is NOT SEG's own published default (12): the first 679-variant run built
# with 12 gave PM4 a real accuracy regression (72->47 recall dropped further at first
# pass, then found to be a false-positive problem specifically) - most PM4 ground-truth
# variants are SHORT in-frame indels (1-4 residues), and padding a 12-residue window out
# from such a short interval mostly measures unrelated flanking sequence, not the
# altered site itself. Ordinary (non-repetitive) human protein sequence's own natural
# amino-acid frequency skew is enough to dip a lucky/unlucky 12-mer under K1=2.2 with no
# real repeat present (confirmed on the actual false positives: GAA/HNF1A/PAH/LDLR/
# CDKL5 windows scored 1.53-2.15 at W=12, all real UniProt data, none of them a genuine
# repeat) - this is exactly the instability SEG's own extend/merge scan exists to filter
# out across a whole-sequence pass, which this single-window use does not run. Widening
# to 24 restores separation on the same real cases (those same windows rise to
# 2.46-2.72 bits at W=24, while RPGR ORF15 and FOXG1's real low-complexity stretches
# stay under 1.1 bits) - an empirical correction for this narrower deployment, not a
# published SEG parameter.
SEG_WINDOW = 24
SEG_TRIGGER_K1 = 2.2
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


def _complexity_window(sequence, start, end):
    """The 1-based [start, end] interval, extended symmetrically to at least
    SEG_WINDOW residues and clamped to the sequence's bounds, read as a substring -
    see SEG_WINDOW's own comment for why this project's value differs from SEG's
    published default. A fixed-length window is scored; an interval shorter than
    that on its own would be scored on too little context to mean anything."""
    length = len(sequence)
    lo, hi = start, end
    grow_low = True
    while hi - lo + 1 < SEG_WINDOW and (lo > 1 or hi < length):
        if grow_low and lo > 1:
            lo -= 1
        elif hi < length:
            hi += 1
        grow_low = not grow_low
    lo, hi = max(1, lo), min(length, hi)
    return sequence[lo - 1:hi]


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
        body, response, used_accession = resolve_covering_entry(self.client, accession, end)
        if body is None:
            return []
        features = body.get("features")
        if not isinstance(features, list):
            features = []
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
        sequence = (body.get("sequence") or {}).get("value")
        complexity = None
        if isinstance(sequence, str) and sequence:
            complexity = sequence_complexity(_complexity_window(sequence, start, end))
            if complexity is not None and complexity <= SEG_TRIGGER_K1:
                repetitive = True
        digest = hashlib.sha256(
            f"{response['body_sha256']}:{start}-{end}".encode("utf-8")).hexdigest()
        return [{
            "category": "region",
            "variant_key": variant.key,
            "evidence_id": f"urn:sha256:{digest}:uniprot-repeat:{transcript}:{start}-{end}",
            "source": self.name,
            "source_version": used_accession,
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
            "sequence_complexity": complexity,
            "isoform_used": used_accession if used_accession != accession else None,
            "policy_note": (
                "UniProt Repeat/Compositional-bias feature overlap, OR Wootton-Federhen/SEG "
                f"sequence complexity <= {SEG_TRIGGER_K1} bits (window >= {SEG_WINDOW} aa), "
                "over the altered protein interval, flagged as an automated first pass "
                "(assessment_method=automated), not a curator's own review. "
                "functional_importance is set only from a Domain/Binding site/Active "
                "site/Motif/Coiled coil (or non-disordered Region) feature overlapping the "
                "same interval, never from sequence complexity alone. isoform_used names an "
                "isoform accession when the canonical entry's own sequence was too short to "
                "cover the altered interval and a longer isoform was used instead - see this "
                "provider's own module docstring and uniprot_features.resolve_covering_entry() "
                "for why, and for why no overlapping repetitive feature or low-complexity "
                "window is read as repetitive=False rather than left unset."
            ),
        }]


__all__ = ["METHOD", "RegionRepeatProvider"]
