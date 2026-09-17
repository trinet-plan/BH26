"""OMIM and Orphanet disease identifiers resolved to MONDO, from MONDO's own mapping set.

[Why this exists]
  PVS1's disease gate compares the disease a case names with the disease a mechanism record
  was curated for, and a comparison of identifier strings only answers when both sides happen
  to use the same vocabulary. A case recorded as OMIM:143890 and a ClinGen curation recorded
  as MONDO:0007254 are the same disease, and without a mapping the gate reads them as
  different and withholds PVS1 for want of evidence that is sitting right there.

[Why only exactMatch]
  MONDO publishes its cross-references as SSSOM with a predicate on every row. Only
  skos:exactMatch asserts the same disease concept; broadMatch and the rest assert a
  relationship, which is what ACMG's disease-match levels keep apart from equivalence on
  purpose. So this reads exactMatch rows and ignores the others rather than treating a
  related concept as the same one.

[Why an ambiguous mapping produces nothing]
  A reverse lookup is only usable if it lands on one MONDO term. In the release this was
  written against, all 10045 OMIM and 9785 Orphanet exactMatch identifiers resolve to exactly
  one MONDO term each, so ambiguity is not a situation that arises today - but a later
  release could introduce one, and picking either term would be choosing a disease on the
  curator's behalf. Such an identifier resolves to nothing and the gate reports it unresolved.

[What it does not do]
  This settles EQUIVALENT only (ACMG disease-match level 2). INCLUDED needs ClinGen's own
  lumping and splitting decisions, and a parent/child relation in the ontology is deliberately
  not equivalence - neither is derivable from this file, and neither is guessed from it.
"""

from __future__ import annotations

import hashlib

MAPPING_URL = "https://purl.obolibrary.org/obo/mondo/mappings/mondo.sssom.tsv"
METHOD = "mondo_sssom_exact_match"
EXACT_MATCH = "skos:exactMatch"
# The vocabularies a case or a curation is realistically recorded in here. Restricting the
# index to these keeps a 13 MB file from becoming a dictionary of everything MONDO has ever
# been cross-referenced to, most of which is not a disease identifier a curator would supply.
RESOLVABLE_PREFIXES = ("OMIM", "Orphanet", "DOID", "MEDGEN")


def parse(text):
    """Return {external id: MONDO id} for unambiguous exactMatch rows, plus the header lines.

    An identifier that resolves to more than one MONDO term is dropped rather than decided.
    """
    lines = text.splitlines()
    header_at = next((index for index, line in enumerate(lines)
                      if line and not line.startswith("#")), None)
    if header_at is None:
        raise ValueError("MONDO mapping set has no header row")
    columns = lines[header_at].split("\t")
    try:
        subject_at = columns.index("subject_id")
        predicate_at = columns.index("predicate_id")
        object_at = columns.index("object_id")
    except ValueError as error:
        raise ValueError("MONDO mapping set is missing a required column") from error

    candidates = {}
    for line in lines[header_at + 1:]:
        fields = line.split("\t")
        if len(fields) <= object_at or fields[predicate_at] != EXACT_MATCH:
            continue
        subject, external = fields[subject_at].strip(), fields[object_at].strip()
        if not subject.startswith("MONDO:") or ":" not in external:
            continue
        if not external.startswith(tuple(f"{prefix}:" for prefix in RESOLVABLE_PREFIXES)):
            continue
        candidates.setdefault(external, set()).add(subject)
    return ({external: next(iter(terms)) for external, terms in candidates.items()
             if len(terms) == 1},
            [line for line in lines[:header_at] if line.startswith("#")])


def version(header_lines):
    """The release this mapping set came from.

    SSSOM carries its provenance in the commented preamble rather than in a dated line of its
    own, so the preamble is hashed: two runs replaying the same file agree, and a file that
    changed underneath a cached result does not silently pass as the same policy.
    """
    if not header_lines:
        return None
    digest = hashlib.sha256("\n".join(header_lines).encode("utf-8")).hexdigest()
    return f"sssom-preamble-{digest[:16]}"


class MondoMappingProvider:
    """Resolves a disease identifier to MONDO, or reports that it could not be resolved."""

    name = "MONDO disease mapping set (SSSOM)"

    def __init__(self, client, *, url=MAPPING_URL):
        self.client = client
        self.url = url

    def normalize(self, condition):
        """A condition mapping for `condition`, or None when it cannot be resolved.

        A MONDO identifier resolves to itself and is reported as such, so a caller can tell a
        case that was already recorded in MONDO from one that had to be mapped - the two are
        different disease-match levels and only the first is an identifier match.
        """
        if not isinstance(condition, str) or ":" not in condition:
            return None
        condition = condition.strip()
        if condition.startswith("MONDO:"):
            return {"input_condition": condition, "normalized_condition": condition,
                    "mapping_type": "identity", "source": self.name,
                    "source_version": None, "retrieved_at": None,
                    "assessment_method": "automated", "method": METHOD}
        response = self.client.fetch(self.url, response_format="text")
        text = response["body"]
        mappings, header = parse(text)
        release = version(header)
        if not release:
            return None
        normalized = mappings.get(condition)
        if not normalized:
            return None
        return {
            "input_condition": condition,
            "normalized_condition": normalized,
            "mapping_type": "equivalent",
            "source": self.name,
            "source_version": release,
            "retrieved_at": response["retrieved_at"],
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": release,
            "response_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        }
