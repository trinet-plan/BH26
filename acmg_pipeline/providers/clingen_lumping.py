"""ClinGen's lumping and splitting decisions: which phenotypes a curated disease covers.

[Why this exists]
  Curations and cases name diseases at different depths, and MONDO ancestry only says the two
  terms are on one path - which is a reason to ask a curator, not an answer. ClinGen's gene
  curation expert panels have already answered it for the diseases they curated: each
  gene-disease validity curation records the MIM phenotypes that were lumped into the disease
  entity and the ones that were deliberately kept out.

  An included phenotype is the curated disease, so a mechanism curated for it applies. An
  excluded one was looked at and ruled out, which is a stronger statement than any
  resemblance: it turns a question a curator would otherwise have to answer into one that has
  already been answered, and it outranks name similarity and ontology ancestry alike.

[Why this reads HTML]
  ClinGen publishes the decision only on the curation page. The gene-validity download names
  the disease, the mode of inheritance and the report URL, but carries no lumping columns;
  there is no JSON representation of the report and no bulk file for it. Only the curations
  for the genes being evaluated are fetched, so a run costs the index plus a page per
  curation rather than the whole corpus.

[Why an unrecognized page yields nothing at all]
  Markup can change without notice, and the dangerous failure is a parse that quietly returns
  an empty included list - a curated disease would then appear to cover nothing and every
  case would be read as excluded from it. So every landmark is required: the report pane, both
  phenotype headings, at least one included phenotype, and the evaluation date. If any is
  absent the page resolves to nothing, and the gate reports the scope as unavailable rather
  than as empty.
"""

from __future__ import annotations

import csv
import html
import io
import re

CURATION_INDEX = "https://search.clinicalgenome.org/kb/gene-validity/download"
METHOD = "clingen_lumping_and_splitting"
REPORT_PANE = 'id="gdvt5"'
INCLUDED_HEADING = "Included MIM Phenotypes"
EXCLUDED_HEADING = "Excluded MIM Phenotypes"
DATE_HEADING = "Evaluation Date"
_MIM = re.compile(r"MIM:(\d{4,8})")
_TAGS = re.compile(r"<[^>]+>")
_SPACE = re.compile(r"\s+")
INDEX_COLUMNS = ("GENE SYMBOL", "DISEASE ID (MONDO)", "MOI", "CLASSIFICATION",
                 "ONLINE REPORT", "CLASSIFICATION DATE")


def parse_index(text):
    """Return one entry per gene-disease validity curation in the published index."""
    rows = list(csv.reader(io.StringIO(text)))
    header_at = next((index for index, row in enumerate(rows)
                      if row and row[0].strip() == "GENE SYMBOL"), None)
    if header_at is None:
        raise ValueError("ClinGen gene-validity index is missing its header row")
    columns = [item.strip() for item in rows[header_at]]
    try:
        at = {name: columns.index(name) for name in INDEX_COLUMNS}
    except ValueError as error:
        raise ValueError("ClinGen gene-validity index is missing a required column") from error
    entries = []
    for row in rows[header_at + 1:]:
        if len(row) <= max(at.values()):
            continue
        values = {name: row[index].strip() for name, index in at.items()}
        gene = values["GENE SYMBOL"]
        # The separator row under the header is all plus signs, not a curation.
        if not gene or gene.startswith("+"):
            continue
        if not values["DISEASE ID (MONDO)"].startswith("MONDO:"):
            continue
        entries.append({
            "gene": gene,
            "condition": values["DISEASE ID (MONDO)"],
            "moi": values["MOI"] or None,
            "classification": values["CLASSIFICATION"] or None,
            "report": values["ONLINE REPORT"] or None,
            "classified_at": values["CLASSIFICATION DATE"] or None,
        })
    return entries


def _text(markup):
    return _SPACE.sub(" ", html.unescape(_TAGS.sub(" ", markup))).strip()


def parse_report(markup):
    """The included and excluded MIM phenotypes on one curation page, or None.

    None means the page did not present the report in the shape this parser recognizes. It is
    never an empty result: a curated disease that appears to include nothing would read as
    excluding every case.
    """
    start = markup.find(REPORT_PANE)
    if start < 0:
        return None
    body = _text(markup[start:])
    positions = [body.find(heading)
                 for heading in (INCLUDED_HEADING, EXCLUDED_HEADING, DATE_HEADING)]
    if any(position < 0 for position in positions) or positions != sorted(positions):
        return None
    included_at, excluded_at, date_at = positions
    included = _MIM.findall(body[included_at:excluded_at])
    excluded = _MIM.findall(body[excluded_at:date_at])
    if not included:
        # Every curation lumps at least the disease it curated, so none of them means the
        # section was not read, not that the disease covers nothing.
        return None
    evaluated = body[date_at + len(DATE_HEADING):date_at + len(DATE_HEADING) + 60]
    evaluated = evaluated.strip(" :|").split(" ")[0].strip()
    if not evaluated:
        return None
    return {
        "included": [f"OMIM:{item}" for item in dict.fromkeys(included)],
        "excluded": [f"OMIM:{item}" for item in dict.fromkeys(excluded)],
        "evaluated_at": evaluated,
    }


class ClinGenLumpingProvider:
    """Builds the disease scope of each curated gene-disease pair, in MONDO identifiers."""

    name = "ClinGen Lumping and Splitting"

    def __init__(self, client, resolver, *, index_url=CURATION_INDEX):
        self.client = client
        # Anything with .normalize(identifier); the MIM phenotypes are published as OMIM.
        self.resolver = resolver
        self.index_url = index_url
        self._index = None

    def index(self):
        if self._index is None:
            response = self.client.fetch(self.index_url, response_format="text")
            self._index = parse_index(response["body"])
        return self._index

    def _mondo(self, identifiers):
        """MONDO terms for the phenotypes that resolve, and the identifiers that did not."""
        resolved, unresolved = [], []
        for item in identifiers:
            mapping = self.resolver.normalize(item)
            if mapping and mapping.get("normalized_condition"):
                resolved.append(mapping["normalized_condition"])
            else:
                unresolved.append(item)
        return sorted(dict.fromkeys(resolved)), unresolved

    def get_scope(self, gene, condition):
        """The disease scope ClinGen recorded for `gene` and `condition`, or None."""
        entry = next((item for item in self.index()
                      if item["gene"] == gene and item["condition"] == condition), None)
        if entry is None or not entry["report"]:
            return None
        response = self.client.fetch(entry["report"], response_format="text")
        report = parse_report(response["body"])
        if report is None:
            return None
        included, included_unresolved = self._mondo(report["included"])
        excluded, excluded_unresolved = self._mondo(report["excluded"])
        return {
            # The curated disease is its own scope, whether or not its MIM entry resolved.
            "included": sorted(dict.fromkeys([condition, *included])),
            "excluded": excluded,
            "source": self.name,
            "source_version": report["evaluated_at"],
            "retrieved_at": response["retrieved_at"],
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": report["evaluated_at"],
            "report": entry["report"],
            "gene_disease_validity": entry["classification"],
            "unresolved_phenotypes": sorted(included_unresolved + excluded_unresolved),
        }
