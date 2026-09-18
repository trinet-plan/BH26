"""TogoVar API population-frequency adapter for normalized GRCh38 variants."""

from decimal import Decimal, InvalidOperation


API_URL = "https://grch38.togovar.org/api/search/variant?stat=0"
API_VERSION = "0.9.1"
MAX_REPORTED_AF_ROUNDING_ERROR = Decimal("0.0005")
FREQUENCY_SOURCE_GROUPS = {
    "gnomad": frozenset({"gnomad_exomes", "gnomad_genomes"}),
    "tommo": frozenset({"tommo"}),
    "jga": frozenset({"jga_snp", "jga_wes", "jga_wgs"}),
    "ncbn": frozenset({"ncbn"}),
    "gem_j": frozenset({"gem_j_wga"}),
    "hgvd": frozenset({"hgvd"}),
}


class TogoVarProvider:
    """Resolve exact GRCh38 alleles through TogoVar without inferring absence."""

    name = "TogoVar"

    def __init__(self, client, *, api_version=API_VERSION, frequency_sources=None):
        if api_version != API_VERSION:
            raise ValueError(f"Only TogoVar API {API_VERSION} is supported")
        requested = list(FREQUENCY_SOURCE_GROUPS) if frequency_sources is None else frequency_sources
        if (not isinstance(requested, list) or not requested
                or not all(isinstance(item, str) for item in requested)
                or len(requested) != len(set(requested))):
            raise ValueError("TogoVar frequency_sources must be a nonempty unique string list")
        unknown = set(requested) - set(FREQUENCY_SOURCE_GROUPS)
        if unknown:
            raise ValueError(f"Unknown TogoVar frequency source groups: {sorted(unknown)}")
        self.client = client
        self.api_version = api_version
        self.frequency_sources = tuple(requested)
        self.datasets = frozenset(
            dataset
            for source in self.frequency_sources
            for dataset in FREQUENCY_SOURCE_GROUPS[source]
        )

    def get_frequency(self, variant, context=None):
        response = self._search_at(variant, variant.pos)
        observations = self._parse_response(response, variant)
        if observations is not None or not self._is_vcf_deletion(variant):
            return observations

        # TogoVar documents a deletion at the first deleted base, with the
        # deleted sequence as REF and no ALT.  A normalized VCF instead
        # anchors the deletion on its preceding retained base.  For example,
        # VCF POS=100 REF=TA ALT=T is TogoVar POS=101 REF=A ALT=<empty>.
        # Do not run this extra lookup for substitutions/insertions: it is
        # only a documented alternate representation of a deletion.
        response = self._search_at(variant, variant.pos + len(variant.alt))
        return self._parse_response(response, variant)

    def _search_at(self, variant, position):
        return self.client.fetch(
            API_URL,
            data={
                "query": {
                    "location": {
                        "chromosome": variant.chrom,
                        "position": position,
                    }
                },
                "limit": 1000,
            },
            dataset_version=f"TogoVar API {self.api_version}",
        )

    def get_frequencies(self, variants):
        """Fetch multiple variants while retaining one replayable cache entry per allele."""
        return {
            variant.key: self.get_frequency(variant)
            for variant in sorted(set(variants), key=lambda item: item.key)
        }

    def _parse_response(self, response, variant):
        body = response.get("body")
        if not isinstance(body, dict) or not isinstance(body.get("data"), list):
            raise ValueError("Unexpected TogoVar response")
        matches = [item for item in body["data"]
                   if isinstance(item, dict) and self._matches_variant(item, variant)]
        if not matches:
            return None
        if len(matches) != 1:
            raise ValueError("TogoVar returned duplicate exact alleles")
        remote = matches[0]
        # The search API normally returns a ``tgv...`` ID, but a valid exact
        # allele can omit it (as observed for RPE65 1-68431328-T-C).  TogoVar
        # documents its coordinate/REF/ALT URL as another canonical variant
        # identifier, so lack of the optional display ID must not discard
        # otherwise complete AC/AN/AF evidence.
        tgv_id = remote.get("id")
        if tgv_id is not None and (not isinstance(tgv_id, str) or not tgv_id.startswith("tgv")):
            raise ValueError("Invalid TogoVar variant ID")
        frequencies = remote.get("frequencies")
        if frequencies is None:
            return None
        if not isinstance(frequencies, list):
            raise ValueError("Invalid TogoVar frequencies")
        observations = [
            self._observation(variant, tgv_id, values, response.get("retrieved_at"))
            for values in frequencies
            if isinstance(values, dict) and values.get("source") in self.datasets
        ]
        return observations or None

    @staticmethod
    def _is_vcf_deletion(variant):
        return len(variant.ref) > len(variant.alt) and variant.ref.startswith(variant.alt)

    @classmethod
    def _matches_variant(cls, item, variant):
        """Match either VCF alleles or TogoVar's documented deletion form."""
        direct = (
            str(item.get("chromosome")), item.get("position"),
            item.get("reference"), item.get("alternate"),
        ) == (variant.chrom, variant.pos, variant.ref, variant.alt)
        if direct:
            return True
        if not cls._is_vcf_deletion(variant):
            return False
        deleted = variant.ref[len(variant.alt):]
        return (
            str(item.get("chromosome")) == variant.chrom
            and item.get("position") == variant.pos + len(variant.alt)
            and item.get("reference") == deleted
            and item.get("alternate") in (None, "")
        )

    def _observation(self, variant, tgv_id, values, retrieved_at):
        if not isinstance(values, dict):
            raise ValueError("Invalid TogoVar frequency observation")
        dataset = values.get("source")
        filters = values.get("filter")
        ac, an, supplied_af = values.get("ac"), values.get("an"), values.get("af")
        if not isinstance(dataset, str) or not dataset:
            raise ValueError("TogoVar frequency source is missing")
        if not isinstance(filters, list) or not all(isinstance(item, str) for item in filters):
            raise ValueError("Invalid TogoVar frequency filters")
        if (not isinstance(ac, int) or isinstance(ac, bool)
                or not isinstance(an, int) or isinstance(an, bool) or an <= 0 or ac < 0 or ac > an):
            raise ValueError("Invalid TogoVar AC/AN")
        try:
            af = Decimal(str(supplied_af))
        except (InvalidOperation, TypeError, ValueError) as exc:
            raise ValueError("Invalid TogoVar AF") from exc
        computed_af = Decimal(ac) / Decimal(an)
        if (not af.is_finite() or not 0 <= af <= 1
                or abs(af - computed_af) > MAX_REPORTED_AF_ROUNDING_ERROR):
            raise ValueError("Inconsistent TogoVar AC/AN/AF")
        quality = "PASS" if filters and set(filters) == {"PASS"} else "FILTERED"
        return {
            "category": "population",
            "variant_key": variant.key,
            "evidence_id": self._evidence_id(variant, tgv_id, dataset),
            "source": self.name,
            "source_version": f"API {self.api_version}",
            "retrieved_at": retrieved_at,
            "population": f"{dataset}:global",
            "AC": ac,
            "AN": an,
            # Some TogoVar datasets expose a rounded AF (for example GEM-J at
            # three decimal places). AC/AN are exact, so downstream thresholds
            # use their quotient and the API value remains available for audit.
            "AF": str(computed_af),
            "reported_AF": str(af),
            "callable": True,
            "quality_status": quality,
            "filters": filters,
            "variant_flags": [],
            "upstream_dataset": dataset,
            "upstream_dataset_version": "NOT_PROVIDED_BY_TOGOVAR_API",
            "togovar_id": tgv_id,
        }

    @staticmethod
    def _evidence_id(variant, tgv_id, dataset):
        """Use TogoVar's exact coordinate URL when its opaque ID is absent."""
        identifier = tgv_id or f"{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}"
        return f"https://grch38.togovar.org/variant/{identifier}#frequency:{dataset}"
