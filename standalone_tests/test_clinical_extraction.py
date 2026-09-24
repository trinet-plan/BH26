from dataclasses import asdict
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from acmg_pipeline.clinical_extraction import extract_clinical_note


note_path = Path("democase/case1_clinical_note_v1.txt")
clinical_note = note_path.read_text(encoding="utf-8")

result = extract_clinical_note(clinical_note)

print(type(result))
print(json.dumps(asdict(result), ensure_ascii=False, indent=2))
