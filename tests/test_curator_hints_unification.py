"""curatorHints unification (2026-09-18, per the user's direction):

- structuredEvidenceItems folded into curatorHints instead of its own
  extension - both are curator-facing disclosures about the same evidence,
  so one list is enough. category="evidence_item" plus extra checked/detail
  keys keep the per-experiment checklist shape a real curator UI doc asks
  for (doc/recs for expert board.docx).
- CuratorHint gained an optional `category` field so a plain hint
  (previously severity+message only) sits in the same shape as the
  automated engine's review/warning/conflict-categorized ones.
"""

import unittest

from acmg_pipeline.common import CuratorHint
from acmg_pipeline.export import _hints_extension, _structured_evidence_item_hints


class Judgment:
    def __init__(self, items):
        self._items = items

    def structured_evidence_items(self):
        return self._items


class StructuredEvidenceItemHintsTests(unittest.TestCase):
    def test_each_item_becomes_an_evidence_item_hint(self):
        judgment = Judgment([{"label": "Yeast complementation: functionally_abnormal",
                               "checked": True, "detail": "growth defect at 37C"}])
        hints = _structured_evidence_item_hints(judgment)
        self.assertEqual(hints, [{
            "severity": "info", "category": "evidence_item",
            "message": "Yeast complementation: functionally_abnormal",
            "checked": True, "detail": "growth defect at 37C",
        }])

    def test_a_judgment_without_the_method_yields_nothing(self):
        self.assertEqual(_structured_evidence_item_hints(object()), [])

    def test_no_items_yields_an_empty_list_not_none(self):
        self.assertEqual(_structured_evidence_item_hints(Judgment([])), [])


class HintsExtensionMergeTests(unittest.TestCase):
    def test_plain_hints_and_evidence_item_hints_share_one_extension(self):
        ext = _hints_extension(
            [CuratorHint("caution", "Only 1 of 2 papers usable")],
            _structured_evidence_item_hints(Judgment([
                {"label": "Assay A", "checked": False, "detail": "inconclusive"},
            ])),
        )
        self.assertEqual(ext.name, "curatorHints")
        self.assertEqual(ext.value, [
            {"severity": "caution", "category": None, "message": "Only 1 of 2 papers usable"},
            {"severity": "info", "category": "evidence_item", "message": "Assay A",
             "checked": False, "detail": "inconclusive"},
        ])

    def test_a_hint_with_an_explicit_category_keeps_it(self):
        ext = _hints_extension([CuratorHint("caution", "msg", category="unevaluated")])
        self.assertEqual(ext.value[0]["category"], "unevaluated")

    def test_nothing_at_all_returns_none(self):
        self.assertIsNone(_hints_extension([], []))


if __name__ == "__main__":
    unittest.main()
