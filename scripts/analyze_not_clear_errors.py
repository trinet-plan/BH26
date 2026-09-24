"""
Error-analysis harness for the literature engine (PS3/BS3/PS4), inspired by
the Evidence Aggregator paper's (Twede et al. 2025, bioRxiv 2025.03.10.642480)
practice of categorizing disagreements by root cause rather than treating
every "not_clear" the same way.

Parses run_integrated_validation_*.py / run_all_demo_variants.py console
logs (the "--- PMID:X ---" / "[Full text] ..." / "LLM's raw judgment:
direction = ..." / "Rationale: ..." block format used throughout this
project) and buckets every not_clear judgment by why the LLM declined,
cross-tabulated against whether full text or only the abstract was
available for that paper.

[First real run, 2026-09-18, 627 per-paper judgments across the 4 demo
cases + all-variant validation logs]
93.6% of judgments were not_clear, and 71.0% of THOSE were "wrong_variant"
(the retrieved paper discusses a different variant) - split almost evenly
between full_text (221) and abstract_only (196) availability, meaning
missing full text is NOT the driver. Most of the remaining
"other_uncategorized" cases were papers that themselves state no
functional validation exists yet, or GWAS/burden studies with no per-
variant data. Conclusion: the high not_clear rate mostly reflects genuine
literature scarcity for often-novel demo variants (search_candidate_pmids
can only return gene-level papers when no variant-specific paper exists),
not a defect in the LLM's judgment - splitting the per-paper LLM call into
a separate "is this variant mentioned" pre-filter would not be expected to
help, since the existing single-call prompt already asks the model to do
that identification first and stop if unsuccessful, and it already does so
correctly in the cases inspected here.
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

BLOCK_RE = re.compile(
    r"--- PMID:(?P<pmid>\d+) ---\n"
    r"\[Full text\] (?P<fulltext_note>.*?)\n"
    r"\[LLM response\].*?\n"
    r"LLM's raw judgment: direction = (?P<direction>\S+)\n"
    r"Rationale: (?P<rationale>.*?)(?=\n\n|\n---|\Z)",
    re.DOTALL,
)

CATEGORY_RULES = [
    ("wrong_variant", re.compile(
        r"different variant|does not (study|report|test|evaluate|mention)|"
        r"not (identified|studied|reported|tested|evaluated|found|mentioned).{0,40}"
        r"(target variant|the variant|this variant)|"
        r"target variant.{0,40}(is not|was not|not)", re.IGNORECASE)),
    ("large_scale_screen_no_specific_value", re.compile(
        r"screen|saturation mutagenesis|deep mutational|VAMP-seq|"
        r"aggregate|does not (state|report).{0,30}(this|the target) variant.{0,20}own", re.IGNORECASE)),
    ("no_cohort_denominator", re.compile(
        r"case report|not a case-control|no comparison group|no denominator|"
        r"single case|not a cohort", re.IGNORECASE)),
    ("wrong_study_type", re.compile(
        r"functional/mechanistic study|knock-in mice|mouse model|"
        r"does not perform.{0,30}functional assay|not applicable for PS3/BS3", re.IGNORECASE)),
    ("conflicting_or_ambiguous", re.compile(
        r"conflict|ambiguous|inconsistent|disagree", re.IGNORECASE)),
]


def categorize(rationale: str) -> str:
    for name, pattern in CATEGORY_RULES:
        if pattern.search(rationale):
            return name
    return "other_uncategorized"


def fulltext_status(note: str) -> str:
    if "Fetched from PubMed" in note and "abstract" not in note.lower():
        return "full_text"
    if "abstract instead" in note.lower():
        return "abstract_only"
    return "unknown"


def analyze(log_paths: list[Path]):
    blocks = []
    for path in log_paths:
        text = path.read_text(errors="replace")
        blocks.extend(m.groupdict() for m in BLOCK_RE.finditer(text))

    not_clear = [b for b in blocks if b["direction"] == "not_clear"]
    print(f"Total per-paper judgments parsed: {len(blocks)}")
    print(f"  not_clear: {len(not_clear)} ({len(not_clear)/len(blocks):.1%})" if blocks else "")

    by_category = Counter()
    by_category_and_fulltext = defaultdict(Counter)
    for b in not_clear:
        cat = categorize(b["rationale"])
        ft = fulltext_status(b["fulltext_note"])
        by_category[cat] += 1
        by_category_and_fulltext[cat][ft] += 1

    print("\n--- not_clear judgments by root-cause category ---")
    for cat, n in by_category.most_common():
        pct = n / len(not_clear) if not_clear else 0
        ft_breakdown = ", ".join(f"{k}={v}" for k, v in by_category_and_fulltext[cat].items())
        print(f"  {cat:36s} {n:4d} ({pct:5.1%})   [{ft_breakdown}]")

    total_ft = Counter(fulltext_status(b["fulltext_note"]) for b in not_clear)
    print("\n--- not_clear judgments by full-text availability (all categories) ---")
    for k, v in total_ft.most_common():
        print(f"  {k:20s} {v:4d} ({v/len(not_clear):.1%})" if not_clear else "")

    print("\n--- 'wrong_variant' rationale samples (first 8) ---")
    wrong_variant = [b for b in not_clear if categorize(b["rationale"]) == "wrong_variant"]
    for b in wrong_variant[:8]:
        print(f"  PMID:{b['pmid']} [{fulltext_status(b['fulltext_note'])}] {b['rationale'][:160]}")


if __name__ == "__main__":
    paths = [Path(p) for p in sys.argv[1:]]
    if not paths:
        print("usage: analyze_not_clear_errors.py <log1> [log2 ...]")
        sys.exit(1)
    analyze(paths)
