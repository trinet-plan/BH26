"""
api_input.py
The top-level Python class for one democase/case*_api_input_*.json file
(e.g. democase/case1_api_input_case1-noise2.json): a raw case-level JSON
envelope with exactly two fields, "vcf" (a synthetic per-variant VCF
string) and "clinical_note" (a free-text clinical narrative).

[Where this sits relative to the eventual API]
  This is upstream of, and distinct from, the "正規化済み構造化データ"
  request-body schema for POST /v1/variants described in
  doc/docker_api_deployment_plan_v1_ja.md section 0-1/2 (still an open TODO
  there, item "このAPIが要求する...入力スキーマを定義する"). These
  api_input.json files are pre-parse demo fixtures, not that normalized
  schema - ApiCaseInput.parse_vcf() (via acmg_pipeline.vcf_record) turns
  the "vcf" string into VariantRecord(s), and each VariantRecord + this
  case's clinical_note is what a future normalization step would map onto
  that eventual per-variant request body.

[Contract: exactly 1 variant per case]
  Per doc section 0-1 ("1リクエストにつき1件POST"), an api_input.json is
  defined to carry exactly 1 variant - `.parse_vcf().record` raises
  ValueError if that's ever violated (e.g. a hand-edited demo file with 2
  data rows). See acmg_pipeline.vcf_record.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from acmg_pipeline.vcf_record import ParsedVcf, parse_vcf


@dataclass
class ApiCaseInput:
    vcf: str
    clinical_note: str

    @staticmethod
    def from_dict(data: dict) -> "ApiCaseInput":
        return ApiCaseInput(vcf=data["vcf"], clinical_note=data.get("clinical_note", ""))

    @staticmethod
    def from_json_file(path: str | Path) -> "ApiCaseInput":
        with open(path, encoding="utf-8") as f:
            return ApiCaseInput.from_dict(json.load(f))

    def parse_vcf(self) -> ParsedVcf:
        return parse_vcf(self.vcf)
