"""dbNSFP prediction scores served by MyVariant.info, pinned to the dbNSFP release.

Ensembl VEP REST returns REVEL and SpliceAI values without saying which release produced
them, and a calibrated threshold is meaningless unless the scored release is known. This
adapter reads the dbNSFP version from the MyVariant.info metadata endpoint and records it on
every score, so a calibration can be matched to the data it was derived from.
"""

import hashlib
from decimal import Decimal, InvalidOperation
from urllib.parse import quote

from acmg_pipeline.providers.http import FetchError, canonical_json


METADATA_URL = "https://myvariant.info/v1/metadata"
VARIANT_URL = "https://myvariant.info/v1/variant"
# Only meta-predictors ClinGen has calibrated recommendations for; SIFT and PolyPhen-2 are
# deliberately not carried forward as criterion input.
PREDICTORS = (
    ("revel", "REVEL", "protein"),
    ("alphamissense", "AlphaMissense", "protein"),
)


def _number(value):
    if isinstance(value, bool) or value is None:
        return None
    try:
        return Decimal(str(value))
    except (InvalidOperation, ValueError):
        return None


def _first(value):
    """dbNSFP reports one entry per transcript, so a field may arrive as a list."""
    if isinstance(value, list):
        values = [item for item in value if item is not None]
        return values[0] if values else None
    return value


class DbnsfpProvider:
    name = "MyVariant.info dbNSFP"

    def __init__(self, client, *, assembly="hg38", release=None):
        if assembly != "hg38":
            raise ValueError("Only GRCh38/hg38 dbNSFP lookups are supported")
        self.client = client
        self.assembly = assembly
        self.release = release
        self.metadata = None

    def dataset_version(self):
        """The dbNSFP release backing the scores, taken from the service metadata."""
        if self.metadata is None:
            response = self.client.fetch(METADATA_URL, dataset_version=self.release)
            body = response["body"]
            source = body.get("src", {}).get("dbnsfp") if isinstance(body, dict) else None
            version = source.get("version") if isinstance(source, dict) else None
            if not version:
                raise ValueError("MyVariant.info metadata does not report a dbNSFP version")
            self.metadata = {
                "dbnsfp_version": version,
                "build_version": body.get("build_version"),
                "retrieved_at": response["retrieved_at"],
                "response_sha256": hashlib.sha256(canonical_json(body).encode()).hexdigest(),
            }
        return self.metadata

    def get_predictions(self, variant, transcript=None):
        """Return versioned, calibration-eligible scores; absence is not a zero score."""
        # dbNSFP scores nonsynonymous SNVs only, and the g. form used here cannot express
        # an indel, so anything else is skipped rather than queried and read as absent.
        if not (len(variant.ref) == 1 and len(variant.alt) == 1
                and {variant.ref, variant.alt} <= set("ACGT")):
            return []
        metadata = self.dataset_version()
        fields = ",".join(f"dbnsfp.{key}" for key, _, _ in PREDICTORS)
        query = quote(f"chr{variant.chrom}:g.{variant.pos}{variant.ref}>{variant.alt}", safe="")
        url = f"{VARIANT_URL}/{query}?assembly={self.assembly}&fields={fields}"
        try:
            response = self.client.fetch(url, dataset_version=metadata["dbnsfp_version"])
        except FetchError as exc:
            if "HTTP_ERROR:404" in str(exc):
                return []
            raise
        body = response["body"]
        dbnsfp = body.get("dbnsfp") if isinstance(body, dict) else None
        if not isinstance(dbnsfp, dict):
            return []
        digest = hashlib.sha256(canonical_json(body).encode()).hexdigest()
        records = []
        for key, predictor, mechanism in PREDICTORS:
            entry = dbnsfp.get(key)
            score = _number(_first(entry.get("score"))) if isinstance(entry, dict) else None
            if score is None:
                continue
            record = {
                "category": "computational", "variant_key": variant.key,
                "evidence_id": f"urn:sha256:{digest}:dbnsfp:{key}",
                "source": self.name,
                "source_version": f'dbNSFP-{metadata["dbnsfp_version"]}',
                "retrieved_at": response["retrieved_at"], "quality_status": "PASS",
                "predictor": predictor,
                "predictor_version": f'dbNSFP-{metadata["dbnsfp_version"]}',
                "mechanism": mechanism, "score": str(score),
                # The scored release is known, so a calibration can be matched to it.
                "calibration_eligible": True,
                "dataset": f'dbNSFP {metadata["dbnsfp_version"]} via MyVariant.info',
                "response_sha256": digest,
                "metadata_sha256": metadata["response_sha256"],
            }
            rankscore = _number(_first(entry.get("rankscore")))
            if rankscore is not None:
                record["rankscore"] = str(rankscore)
            prediction = _first(entry.get("pred"))
            if isinstance(prediction, str):
                record["classification"] = prediction
            if transcript:
                record["transcript"] = transcript
            records.append(record)
        return records
