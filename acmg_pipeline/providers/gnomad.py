"""gnomAD GraphQL population-frequency adapter for normalized GRCh38 variants."""

from decimal import Decimal


QUERY = """
query VariantFrequency($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id reference_genome chrom pos ref alt flags
    exome { ac an af homozygote_count hemizygote_count filters populations { id ac an homozygote_count hemizygote_count } }
    genome { ac an af homozygote_count hemizygote_count filters populations { id ac an homozygote_count hemizygote_count } }
  }
}
""".strip()

PRIMARY_POPULATIONS = {"afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas"}


class GnomadProvider:
    name = "gnomAD"

    def __init__(self, client, *, dataset="gnomad_r4", release="4.1.1"):
        if dataset != "gnomad_r4" or not release:
            raise ValueError("Only version-recorded gnomAD r4 GRCh38 data are supported")
        self.client = client
        self.dataset = dataset
        self.release = release

    def get_frequency(self, variant, context=None):
        variant_id = f"{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}"
        response = self.client.fetch(
            "https://gnomad.broadinstitute.org/api",
            data={"query": QUERY, "variables": {"variantId": variant_id,
                                                  "dataset": self.dataset}},
            dataset_version=self.release, allow_application_errors=True,
        )
        return self._parse_response(response, variant, "variant")

    def get_frequencies(self, variants):
        """Fetch a deterministic batch without treating absent variants as zero-frequency."""
        variants = sorted(set(variants), key=lambda item: item.key)
        if not variants:
            return {}
        declarations = ["$dataset: DatasetId!"]
        fields = []
        variables = {"dataset": self.dataset}
        selection = QUERY.split("variant(variantId:", 1)[1].split("{", 1)[1].rsplit("}", 2)[0]
        for index, variant in enumerate(variants):
            alias = f"v{index}"
            variable = f"variantId{index}"
            declarations.append(f"${variable}: String!")
            variables[variable] = self._variant_id(variant)
            fields.append(
                f"{alias}: variant(variantId: ${variable}, dataset: $dataset) {{{selection}}}"
            )
        query = f"query VariantFrequencies({', '.join(declarations)}) {{ {' '.join(fields)} }}"
        response = self.client.fetch(
            "https://gnomad.broadinstitute.org/api",
            data={"query": query, "variables": variables},
            dataset_version=self.release, allow_application_errors=True,
        )
        return {
            variant.key: self._parse_response(response, variant, f"v{index}")
            for index, variant in enumerate(variants)
        }

    def _parse_response(self, response, variant, field):
        body = response["body"]
        if not isinstance(body, dict):
            raise ValueError("Unexpected gnomAD response")
        data = body.get("data")
        remote = data.get(field) if isinstance(data, dict) else None
        all_errors = [item for item in body.get("errors", []) if isinstance(item, dict)]
        field_errors = [item for item in all_errors if item.get("path") == [field]]
        unscoped_errors = [item for item in all_errors if item.get("path") is None]
        fatal_unscoped = [item for item in unscoped_errors
                          if "not found" not in item.get("message", "").lower()]
        if fatal_unscoped:
            raise ValueError("gnomAD query failed")
        if remote is None:
            messages = [item.get("message", "") for item in field_errors or unscoped_errors]
            if messages and all("not found" in message.lower() for message in messages):
                return None
            raise ValueError("gnomAD query failed")
        if field_errors:
            raise ValueError("Partial gnomAD response")
        variant_id = self._variant_id(variant)
        returned = (str(remote.get("chrom")), remote.get("pos"), remote.get("ref"),
                    remote.get("alt"), remote.get("reference_genome"))
        if returned != (variant.chrom, variant.pos, variant.ref, variant.alt, "GRCh38"):
            raise ValueError("gnomAD response variant mismatch")
        if remote.get("variant_id") != variant_id:
            raise ValueError("gnomAD response ID mismatch")
        output = []
        variant_flags = remote.get("flags") or []
        if not isinstance(variant_flags, list):
            raise ValueError("Invalid gnomAD flags")
        for callset in ("exome", "genome"):
            values = remote.get(callset)
            if values is None:
                continue
            if not isinstance(values, dict):
                raise ValueError("Invalid gnomAD callset")
            filters = values.get("filters") or []
            if not isinstance(filters, list):
                raise ValueError("Invalid gnomAD filters")
            output.append(self._observation(variant, callset, "global", values,
                                            filters, variant_flags, response["retrieved_at"]))
            populations = values.get("populations") or []
            if not isinstance(populations, list):
                raise ValueError("Invalid gnomAD populations")
            for population in populations:
                if not isinstance(population, dict) or not population.get("id"):
                    raise ValueError("Invalid gnomAD population observation")
                if population["id"] not in PRIMARY_POPULATIONS:
                    continue
                output.append(self._observation(variant, callset, population["id"], population,
                                                filters, variant_flags, response["retrieved_at"]))
        return output or None

    @staticmethod
    def _variant_id(variant):
        return f"{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}"

    def _observation(self, variant, callset, population, values, filters, flags, retrieved_at):
        ac, an = values.get("ac"), values.get("an")
        if not isinstance(ac, int) or isinstance(ac, bool) or not isinstance(an, int) or an <= 0:
            raise ValueError("Invalid gnomAD AC/AN")
        computed_af = Decimal(ac) / Decimal(an)
        supplied_af = values.get("af")
        if supplied_af is not None and abs(Decimal(str(supplied_af)) - computed_af) > Decimal("0.000001"):
            raise ValueError("Inconsistent gnomAD AC/AN/AF")
        quality = "PASS" if not filters and not flags else "FILTERED"
        variant_id = f"{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}"
        # homozygote_count/hemizygote_count (BS2's own signal - added 2026-09-25, unused by
        # BA1/BS1/PM2) are read the same way AC/AN are: an int when gnomAD reports one, and
        # left absent (None) rather than coerced to 0 when the field is missing, so BS2 can
        # tell "this source has no genotype-count data" from "this source counted zero".
        homozygote_count = values.get("homozygote_count")
        hemizygote_count = values.get("hemizygote_count")
        for label, value in (("homozygote_count", homozygote_count), ("hemizygote_count", hemizygote_count)):
            if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 0):
                raise ValueError(f"Invalid gnomAD {label}")
        return {
            "category": "population", "variant_key": variant.key,
            "evidence_id": (f"https://gnomad.broadinstitute.org/variant/{variant_id}"
                            f"?dataset={self.dataset}#{callset}:{population}"),
            "source": self.name, "source_version": self.release,
            "retrieved_at": retrieved_at, "population": f"{callset}:{population}",
            "AC": ac, "AN": an, "AF": str(computed_af), "callable": True,
            "quality_status": quality, "filters": filters, "variant_flags": flags,
            "homozygote_count": homozygote_count, "hemizygote_count": hemizygote_count,
        }
