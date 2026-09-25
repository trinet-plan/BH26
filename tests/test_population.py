from acmg_pipeline.constants import CriterionStatus
import unittest
from types import SimpleNamespace

from acmg_pipeline.automated_core.interface import inputs_from_prepared_record
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria import ba1, bs1, bs2, pm2
from acmg_pipeline.providers.local import LocalPopulationProvider
from acmg_pipeline.services.evidence import EvidenceService
from acmg_pipeline.services.population import PopulationService
from acmg_pipeline.automated_va_spec import to_evidence_line, validate_1_0_1


class PopulationTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.input = {"variant": self.variant.to_dict()}
        self.config = {code: {"minimum_an": 2000, "policy_source": "synthetic-test-policy",
                              "policy_version": "1"} for code in ("PM2", "BA1", "BS1", "BS2")}
        # "Absent from controls" is stated, not defaulted: pm2.evaluate() reports an unset
        # max_af as an unconfigured policy rather than running as the strictest threshold.
        self.config["PM2"]["max_af"] = 0
        # 5% is the ACMG/AMP 2015 general figure and a point-estimate cutoff, so it is
        # compared with AF; a VCEP specification replaces the pair, not just the number.
        self.config["BA1"].update({"max_af": 0.05, "frequency_statistic": "af"})
        # BS1 falls back to this when no disease-specific threshold applies. Like every other
        # threshold in this project it lives in the policy, not in the code.
        self.config["BS1"].update({"default_max_credible_af": 0.0001,
                                   "default_source": "synthetic-default-policy",
                                   "default_source_version": "1",
                                   "default_frequency_statistic": "faf95",
                                   "default_comparison": ">"})
        # BS2's configured default max_count, the same "used unless a gene-specific
        # override replaces it entirely" shape as BA1's max_af.
        self.config["BS2"]["max_count"] = 0

    def observation(self, ac=0, an=10000, **extra):
        return {"variant_key": self.variant.key, "evidence_id": "test:frequency",
                "source": "synthetic", "source_version": "1", "population": "TEST",
                "retrieved_at": "2026-09-14T00:00:00Z", "AC": ac, "AN": an,
                "AF": ac / an, "quality_status": "PASS", "callable": True,
                "homozygote_count": None, "hemizygote_count": None, **extra}

    def services(self, *batches, evidence_records=()):
        return SimpleNamespace(
            population=PopulationService([
                LocalPopulationProvider(f"test-{i}", batch) for i, batch in enumerate(batches)]),
            evidence=EvidenceService(list(evidence_records)),
        )

    def evaluate(self, module, services):
        variant, clinical_note = inputs_from_prepared_record(self.input)
        return module.evaluate(variant, clinical_note, services, self.config)

    def test_absence_is_supporting(self):
        value = self.evaluate(pm2, self.services([self.observation()]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evidence_outcome, "PM2_supporting")

    def test_no_record_is_a_flagged_prediction_not_af_zero(self):
        """No record at all is not read as a confirmed AF 0, but (2026-09-17 policy change)
        reaches a flagged supporting-strength MET rather than withholding a call."""
        value = self.evaluate(pm2, self.services([]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.strength, "supporting")

    def test_unset_threshold_is_not_the_strictest_threshold(self):
        """An absent max_af used to default to 0 and score PM2 on an unstated policy."""
        del self.config["PM2"]["max_af"]
        value = self.evaluate(pm2, self.services([self.observation()]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("PM2.max_af", value.missing_inputs)

    def test_configured_zero_threshold_still_scores(self):
        """max_af = 0 is the standard "absent from controls" policy, not a missing one."""
        self.config["PM2"]["max_af"] = 0
        value = self.evaluate(pm2, self.services([self.observation()]))
        self.assertEqual(value.status, CriterionStatus.MET)

    def test_high_regional_af_is_not_hidden(self):
        value = self.evaluate(pm2, self.services(
            [self.observation()], [self.observation(ac=100, population="REGIONAL")]
        ))
        self.assertEqual(value.status, CriterionStatus.NOT_MET)

    def test_unavailable_provider_prevents_rarity(self):
        value = self.evaluate(pm2, self.services([self.observation()], []))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)

    def test_excluded_small_subgroup_does_not_invalidate_usable_population(self):
        usable = self.observation()
        too_small = self.observation(an=100, population="SMALL")
        self.assertEqual(self.evaluate(pm2, self.services([usable, too_small])).status, CriterionStatus.MET)
        self.assertEqual(self.evaluate(ba1, self.services([usable, too_small])).status, CriterionStatus.NOT_MET)

    def test_low_an_and_unknown_callability(self):
        for item in (self.observation(an=100), self.observation(callable=False),
                     self.observation(AF="NaN"), self.observation(AF=0.3)):
            with self.subTest(item=item):
                value = self.evaluate(pm2, self.services([item]))
                self.assertEqual(value.status, CriterionStatus.UNKNOWN)

    def test_ba1_boundary_and_exception(self):
        self.input["ba1_exception_assessment"] = {
            "source": "test", "source_version": "1", "reviewed_at": "2026-09-14",
            "is_exception": False,
        }
        for ac, expected in ((499, CriterionStatus.NOT_MET), (500, CriterionStatus.NOT_MET), (501, CriterionStatus.MET)):
            value = self.evaluate(ba1, self.services([self.observation(ac)]))
            self.assertEqual(value.status, expected)
        self.input["ba1_exception_assessment"]["is_exception"] = True
        self.assertEqual(
            self.evaluate(ba1, self.services([self.observation(501)])).status,
            CriterionStatus.NOT_MET,
        )

    def test_ba1_reports_an_unconfigured_threshold_instead_of_assuming_five_percent(self):
        """5% was hardcoded, so a run silently disagreed with any VCEP that replaces it."""
        for key in ("max_af", "frequency_statistic"):
            del self.config["BA1"][key]
        value = self.evaluate(ba1, self.services([self.observation(600)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["BA1.max_af", "BA1.frequency_statistic"])

    def test_ba1_honours_a_vcep_threshold_below_the_general_figure(self):
        """The ClinGen Cardiomyopathy VCEP uses 0.1%, not the ACMG/AMP general 5%."""
        self.config["BA1"]["max_af"] = 0.001
        self.input["ba1_exception_assessment"] = {
            "source": "test", "source_version": "1", "reviewed_at": "2026-09-14",
            "is_exception": False}
        value = self.evaluate(ba1, self.services([self.observation(100)]))  # AF 1%
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["stand_alone_threshold"], "0.001")

    def test_ba1_uses_the_statistic_the_policy_names(self):
        observation = self.observation(516, an=10000)  # AF 0.0516, FAF95 0.0480
        self.input["ba1_exception_assessment"] = {
            "source": "test", "source_version": "1", "reviewed_at": "2026-09-14",
            "is_exception": False}
        for statistic, expected in (("af", CriterionStatus.MET), ("faf95", CriterionStatus.NOT_MET)):
            with self.subTest(statistic=statistic):
                self.config["BA1"]["frequency_statistic"] = statistic
                value = self.evaluate(ba1, self.services([observation]))
                self.assertEqual(value.status, expected)
                self.assertEqual(value.provenance["frequency_statistic"], statistic)

    def test_ba1_rejects_a_threshold_outside_zero_to_one(self):
        for bad in (5, 0, "high"):
            with self.subTest(max_af=bad):
                self.config["BA1"]["max_af"] = bad
                value = self.evaluate(ba1, self.services([self.observation(600)]))
                self.assertEqual(value.status, CriterionStatus.UNKNOWN)
                self.assertEqual(value.missing_inputs, ["BA1.max_af"])

    def test_ba1_stops_on_an_unconfigured_population_policy(self):
        """A stand-alone threshold does not substitute for the population quality policy."""
        del self.config["BA1"]["minimum_an"]
        value = self.evaluate(ba1, self.services([self.observation(600)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["BA1.minimum_an"])

    def test_ba1_rejects_a_statistic_it_cannot_compute(self):
        self.config["BA1"]["frequency_statistic"] = "popmax"
        value = self.evaluate(ba1, self.services([self.observation(600)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["BA1.frequency_statistic"])

    def test_ba1_exception_list_is_policy_not_an_evidence_item(self):
        self.input["ba1_exception_assessment"] = {
            "source": "ClinGen SVI BA1 exception list", "source_version": "2018",
            "reviewed_at": "2026-09-15", "is_exception": False}
        value = self.evaluate(ba1, self.services([self.observation(501)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertTrue(all(item.get("evidence_id") for item in value.evidence))
        self.assertFalse(value.provenance["ba1_exception_assessment"]["is_exception"])
        line = to_evidence_line(value)
        validate_1_0_1(line, "BA1")
        self.assertEqual(line["directionOfEvidenceProvided"], "disputes")

    def test_ba1_gene_specific_override_replaces_the_configured_default(self):
        """automated_core.context's gene_frequency_thresholds feeds
        input_data["ba1_threshold_override"], which criteria/ba1.py's threshold_policy()
        checks before config["BA1"] - a real ClinGen RAF1 RASopathy VCEP-style number
        (0.05%) replaces the configured 5% default entirely, not just when it is stricter."""
        self.input["ba1_threshold_override"] = {
            "max_af": 0.0005, "frequency_statistic": "faf95", "comparison": ">=",
            "source": "ClinGen RASopathy VCEP RAF1", "source_version": "2.3",
            "reviewed_at": "2026-09-22"}
        self.input["ba1_exception_assessment"] = {
            "source": "test", "source_version": "1", "reviewed_at": "2026-09-14",
            "is_exception": False}
        value = self.evaluate(ba1, self.services([self.observation(60, an=10000)]))  # AF 0.6%
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["threshold_scope"], "gene_specific")
        self.assertEqual(value.provenance["stand_alone_threshold"], "0.0005")
        self.assertEqual(value.provenance["comparison"], ">=")

    def test_ba1_honours_the_gene_specific_comparison_operator(self):
        """">" and ">=" are both real ClinGen VCEP wordings (see services/population.py's own
        COMPARISONS docstring) - a variant sitting exactly on the boundary must be MET under
        ">=" and NOT_MET under strict ">"."""
        self.input["ba1_exception_assessment"] = {
            "source": "test", "source_version": "1", "reviewed_at": "2026-09-14",
            "is_exception": False}
        for symbol, expected in ((">=", CriterionStatus.MET), (">", CriterionStatus.NOT_MET)):
            with self.subTest(symbol=symbol):
                self.input["ba1_threshold_override"] = {
                    "max_af": 0.001, "frequency_statistic": "af", "comparison": symbol,
                    "source": "test", "source_version": "1", "reviewed_at": "2026-09-22"}
                value = self.evaluate(ba1, self.services([self.observation(10, an=10000)]))  # AF exactly 0.001
                self.assertEqual(value.status, expected)

    def test_an_incomplete_gene_specific_override_falls_back_to_the_default(self):
        """A partial override (missing provenance) is not applied - config["BA1"]'s own
        complete default is used instead, same as if no override existed at all."""
        self.input["ba1_threshold_override"] = {"max_af": 0.0005, "frequency_statistic": "faf95"}
        self.input["ba1_exception_assessment"] = {
            "source": "test", "source_version": "1", "reviewed_at": "2026-09-14",
            "is_exception": False}
        value = self.evaluate(ba1, self.services([self.observation(600)]))
        self.assertEqual(value.provenance["threshold_scope"], "default")

    def test_bs1_requires_matching_disease(self):
        """A matching threshold is used as itself, and is reported as disease-specific."""
        self.input.update({"condition": "test:disease", "inheritance": "autosomal_dominant",
                           "disease_frequency_threshold": {
                               "condition": "test:disease", "inheritance": "autosomal_dominant",
                               "max_credible_af": 0.001, "source": "test", "source_version": "1",
                               "reviewed_at": "2026-09-14"}})
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["disease_frequency_threshold"]["max_credible_af"], 0.001)
        self.assertEqual(value.provenance["threshold_scope"], "disease_specific")
        # The curated threshold does not name the statistic it was calibrated against, so it
        # is read as a point-estimate cutoff and the assumption is surfaced, not buried.
        self.assertEqual(value.provenance["frequency_statistic"], "af")
        # Neither the statistic nor the comparison is named, so both are read the weaker way
        # and both assumptions are surfaced.
        self.assertEqual(value.provenance["comparison"], ">")
        self.assertEqual(value.review_points, [
            "Confirm the curated threshold is a point-estimate AF cutoff, "
            "not a filtering allele frequency",
            "Confirm the specification says allele frequency > the threshold, not >=",
        ])

    def test_bs1_uses_the_statistic_the_threshold_names(self):
        """AF and FAF are different numbers, so a threshold is compared with its own one."""
        observation = self.observation(12, an=15256)  # AF 0.000787, FAF95 0.000491
        for statistic, expected in (("af", CriterionStatus.MET), ("faf95", CriterionStatus.NOT_MET)):
            with self.subTest(statistic=statistic):
                self.with_threshold(max_credible_af=0.0006, frequency_statistic=statistic)
                value = self.evaluate(bs1, self.services([observation]))
                self.assertEqual(value.status, expected)
                self.assertEqual(value.provenance["frequency_statistic"], statistic)

    def test_bs1_applies_the_comparison_the_specification_writes(self):
        """A variant sitting exactly on the threshold is the one the panel drew a line at.

        ClinGen's PAH and Lysosomal Storage Disorders panels write BS1 as AF > their
        threshold; Hearing Loss, PTEN, Glaucoma and Familial Hypercholesterolemia write >=.
        """
        observation = self.observation(10, an=10000)  # AF exactly 0.001
        for symbol, expected in ((">", CriterionStatus.NOT_MET), (">=", CriterionStatus.MET)):
            with self.subTest(comparison=symbol):
                self.with_threshold(max_credible_af=0.001, frequency_statistic="af",
                                    comparison=symbol)
                value = self.evaluate(bs1, self.services([observation]))
                self.assertEqual(value.status, expected)
                self.assertEqual(value.provenance["comparison"], symbol)

    def test_bs1_assumes_the_weaker_comparison_and_says_so(self):
        """Assuming ">" can only withhold BS1, never assert it on a boundary variant."""
        self.with_threshold(max_credible_af=0.001, frequency_statistic="af", comparison=None)
        value = self.evaluate(bs1, self.services([self.observation(10, an=10000)]))
        self.assertEqual(value.status, CriterionStatus.NOT_MET)
        self.assertEqual(value.provenance["comparison"], ">")
        self.assertIn("Confirm the specification says allele frequency > the threshold, not >=",
                      value.review_points)

    def test_bs1_rejects_a_comparison_it_cannot_apply(self):
        self.with_threshold(comparison="approximately")
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["disease_frequency_threshold.comparison"])

    def test_bs1_rejects_a_statistic_it_cannot_compute(self):
        self.with_threshold(frequency_statistic="popmax")
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs,
                         ["disease_frequency_threshold.frequency_statistic"])

    # BS1's maximum credible frequency is a per-disease policy, so it is reached through a
    # chain of gates before any observation is compared. Each gate is a different job for a
    # curator - supply a disease, retarget a threshold, record its provenance, fix a bad
    # value - and reporting them apart is the whole point of the chain, so each is pinned.
    THRESHOLD = {"condition": "test:disease", "inheritance": "autosomal_dominant",
                 "max_credible_af": 0.001, "source": "test", "source_version": "1",
                 "reviewed_at": "2026-09-14"}

    def with_threshold(self, **overrides):
        """Prepare a record whose curated threshold differs from THRESHOLD in one way."""
        threshold = {**self.THRESHOLD, **overrides}
        for key, value in list(threshold.items()):
            if value is None:
                del threshold[key]
        self.input.update({"condition": "test:disease", "inheritance": "autosomal_dominant",
                           "disease_frequency_threshold": threshold})

    def assertFellBackToDefault(self, value, *, because):
        """A default-threshold verdict must say so, and must not pass as disease-specific."""
        self.assertEqual(value.provenance["threshold_scope"], "default")
        self.assertIn(because, value.summary)
        self.assertIn("not a disease-specific threshold", value.summary)
        self.assertIn("Confirm BS1 against a disease-specific maximum credible frequency",
                      value.review_points)

    def test_bs1_without_a_disease_context_uses_the_default_threshold(self):
        """The default threshold is not approved for any disease, so an exceedance under it
        reaches MET only as a flagged prediction (2026-09-17 policy change: unapproved no
        longer withholds the call, it travels in review_points/provenance instead)."""
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertFellBackToDefault(value, because="no disease context was supplied")

    def test_bs1_without_a_curated_threshold_uses_the_default_threshold(self):
        self.input["condition"] = "test:disease"
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertFellBackToDefault(value, because="no curated maximum credible allele frequency")

    def test_bs1_threshold_for_another_disease_does_not_apply_to_this_one(self):
        self.with_threshold(condition="test:other-disease")
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertFellBackToDefault(value, because="test:other-disease")

    def test_bs1_reports_an_unconfigured_default_instead_of_inventing_one(self):
        """No disease-specific threshold and no configured default is not a verdict."""
        for key in ("default_max_credible_af", "default_source", "default_source_version",
                    "default_frequency_statistic"):
            del self.config["BS1"][key]
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["BS1.default_max_credible_af",
                                                "BS1.default_source",
                                                "BS1.default_source_version",
                                                "BS1.default_frequency_statistic"])

    def test_bs1_rejects_a_default_outside_zero_to_one(self):
        self.config["BS1"]["default_max_credible_af"] = 5
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["BS1.default_max_credible_af"])

    def test_bs1_rejects_a_threshold_outside_zero_to_one(self):
        for bad in (5, 0, -0.1, "high"):
            with self.subTest(max_credible_af=bad):
                self.with_threshold(max_credible_af=bad)
                value = self.evaluate(bs1, self.services([self.observation(100)]))
                self.assertEqual(value.status, CriterionStatus.UNKNOWN)
                self.assertEqual(value.missing_inputs,
                                 ["disease_frequency_threshold.max_credible_af"])

    def test_bs1_names_only_the_provenance_fields_actually_absent(self):
        """An unprovenanced threshold is not disease-specific, and the summary says which
        fields would make it so - `inheritance`, which is present, is not named."""
        self.with_threshold(source=None, reviewed_at=None)
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertFellBackToDefault(value, because="does not record source, reviewed_at")
        self.assertNotIn("inheritance", value.summary)

    def test_bs1_will_not_apply_a_threshold_from_another_inheritance_mode(self):
        self.with_threshold(inheritance="autosomal_recessive")
        self.input["inheritance"] = "autosomal_dominant"
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertFellBackToDefault(value, because="'autosomal_recessive' inheritance")

    def test_bs1_stops_on_an_unconfigured_population_policy(self):
        """A disease threshold does not substitute for the population quality policy."""
        self.with_threshold()
        del self.config["BS1"]["minimum_an"]
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["BS1.minimum_an"])

    def test_bs1_below_threshold_with_every_source_resolved_is_not_met(self):
        self.with_threshold()
        value = self.evaluate(bs1, self.services([self.observation(1)]))
        self.assertEqual(value.status, CriterionStatus.NOT_MET)
        self.assertIsNone(value.strength)

    def test_bs1_below_threshold_with_an_unreachable_source_stays_unknown(self):
        """An incomplete search is not evidence that no population exceeds the threshold."""
        self.with_threshold()
        # A second provider that holds nothing for this variant is recorded as a failure
        # (NO_OBSERVATION), which is what an unreachable source looks like to BS1.
        value = self.evaluate(bs1, self.services([self.observation(1)], []))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("complete_population_evidence", value.missing_inputs)

    def test_bs1_threshold_is_policy_not_an_evidence_item(self):
        """Every exported evidence item needs an IRI; the threshold has none, so it stays out."""
        self.input.update({"condition": "test:disease", "inheritance": "autosomal_dominant",
                           "disease_frequency_threshold": {
                               "condition": "test:disease", "inheritance": "autosomal_dominant",
                               "max_credible_af": 0.001, "source": "test", "source_version": "1",
                               "reviewed_at": "2026-09-14"}})
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertTrue(value.evidence)
        self.assertTrue(all(item.get("evidence_id") for item in value.evidence))
        line = to_evidence_line(value)
        validate_1_0_1(line, "BS1")
        self.assertEqual(line["evidenceOutcome"]["primaryCoding"]["code"], "BS1")
        self.assertEqual(line["directionOfEvidenceProvided"], "disputes")

    def test_bs2_unsupported_inheritance_is_unknown(self):
        """AD/XLD/mitochondrial are a deliberate scope boundary (see bs2.py's own module
        docstring): a heterozygous observation would double-count BA1/BS1's own signal."""
        self.input["inheritance"] = "autosomal_dominant"
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=5)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("autosomal_dominant", value.summary)

    def test_bs2_no_inheritance_is_unknown(self):
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=5)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["inheritance"])

    def test_bs2_semidominant_uses_homozygote_count(self):
        """Real ClinGen Gene-Disease Validity data (2026-09-25): LDLR/familial
        hypercholesterolemia is curated as "SD" (semidominant), not AD or AR - heterozygotes
        have a milder/later-onset phenotype, homozygotes the severe early-onset one BS2's own
        question is actually about, so this reads the same field as autosomal_recessive."""
        self.input["inheritance"] = "semidominant"
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=1)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["zygosity_field"], "homozygote_count")

    def test_bs2_semidominant_abbreviation_normalizes(self):
        """"SD", ClinGen Gene-Disease Validity's own real MOI abbreviation, must resolve the
        same way as the spelled-out form."""
        self.input["inheritance"] = "SD"
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=1)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["inheritance"], "semidominant")

    def gene_disease(self, inheritance=None, moi=None, condition=None):
        """One Gene2Phenotype ("gene_disease") record and one ClinGen Gene-Disease Validity
        ("gene_disease_validity") record for self.variant, each carrying its own inheritance
        field under its own real name (inheritance vs moi) - the two categories bs2.py's
        _resolve_inheritance() reads, mirroring the shape pvs1.py's own _candidate_
        conditions() already consumes from these same two categories."""
        records = []
        if inheritance is not None:
            records.append({
                "category": "gene_disease", "variant_key": self.variant.key,
                "evidence_id": "test:gene-disease", "source": "test", "source_version": "1",
                "retrieved_at": "2026-09-25T00:00:00Z", "quality_status": "PASS",
                "assessment_method": "automated", "method": "test", "policy_version": "1",
                "inheritance": inheritance, "condition": condition,
            })
        if moi is not None:
            records.append({
                "category": "gene_disease_validity", "variant_key": self.variant.key,
                "evidence_id": "test:gene-disease-validity", "source": "test",
                "source_version": "1", "retrieved_at": "2026-09-25T00:00:00Z",
                "quality_status": "PASS", "moi": moi, "condition": condition,
            })
        return records

    def test_bs2_resolves_inheritance_from_agreeing_gene_disease_evidence(self):
        """No input_data["inheritance"] at all (the real ERepo-sourced-record gap this was
        added for, 2026-09-25) - both curated sources agree, so BS2 can still run."""
        value = self.evaluate(bs2, self.services(
            [self.observation(homozygote_count=1)],
            evidence_records=self.gene_disease(inheritance="AR", moi="Autosomal recessive")))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["inheritance"], "autosomal_recessive")
        self.assertEqual(value.provenance["inheritance_source"],
                         "resolved_from_gene_disease_evidence")

    def test_bs2_does_not_guess_when_gene_disease_sources_disagree(self):
        """A RYR1-style gene curated for more than one disease with different inheritance
        modes must not have one auto-picked for it - see _resolve_inheritance()'s own
        docstring for why this mirrors pvs1.py's "offered, never guessed" caution."""
        value = self.evaluate(bs2, self.services(
            [self.observation(homozygote_count=5)],
            evidence_records=self.gene_disease(inheritance="AD", moi="Autosomal recessive")))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["inheritance"])

    def test_bs2_a_supplied_inheritance_is_not_overridden_by_resolution(self):
        """input_data["inheritance"], when present, is used as-is - _resolve_inheritance()
        is only a fallback, never a second opinion overriding a value already supplied."""
        self.input["inheritance"] = "autosomal_recessive"
        value = self.evaluate(bs2, self.services(
            [self.observation(homozygote_count=1)],
            evidence_records=self.gene_disease(inheritance="AD", moi="AD")))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["inheritance_source"], "supplied")

    def test_bs2_agreement_across_different_mondo_ids_for_the_same_gene_still_resolves(self):
        """Real ERepo data (2026-09-25): PDHA1 carries 6 real curated records unanimously
        naming x_linked inheritance, spanning 4 DIFFERENT MONDO IDs for what every source
        agrees is the same underlying disease at different granularity - exact condition
        matching used to throw this unanimous agreement away entirely (see
        _resolve_inheritance()'s own docstring for the reverted approach and why). Uses
        autosomal_recessive/homozygote_count here (not PDHA1's own x_linked) only to reuse
        this class's default autosomal test variant/observation."""
        self.input["condition"] = "MONDO:0000009"
        value = self.evaluate(bs2, self.services(
            [self.observation(homozygote_count=1)],
            evidence_records=(
                self.gene_disease(inheritance="AR", moi="Autosomal recessive",
                                  condition="MONDO:0000001")
                + self.gene_disease(inheritance="AR", moi="Autosomal recessive",
                                    condition="MONDO:0000002"))))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["inheritance"], "autosomal_recessive")

    def test_bs2_recessive_homozygote_above_default_threshold_is_a_draft_met(self):
        """The only threshold this project has today is config["BS2"]'s configured default
        (no curated bs2_threshold_override entries exist yet), so every MET is currently a
        flagged draft prediction, mirroring bs1.py's own default-scope MET-while-DRAFT
        pattern."""
        self.input["inheritance"] = "autosomal_recessive"
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=1)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.strength, "strong")
        self.assertEqual(value.provenance["policy_status"], "DRAFT")
        self.assertEqual(value.provenance["zygosity_field"], "homozygote_count")
        self.assertIn("prediction pending curator sign-off", value.summary)
        self.assertIn("Approve or replace the draft BS2 genotype-count threshold before "
                      "BS2 is used in a classification", value.review_points)

    def test_bs2_recessive_at_or_below_threshold_is_not_met(self):
        self.input["inheritance"] = "autosomal_recessive"
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=0)]))
        self.assertEqual(value.status, CriterionStatus.NOT_MET)

    def test_bs2_x_linked_uses_hemizygote_count_on_chrx(self):
        variant = Variant("GRCh38", "X", 2, "C", "T")
        self.input = {"variant": variant.to_dict(), "inheritance": "x_linked_recessive"}

        def observation(ac=0, an=10000, **extra):
            return {"variant_key": variant.key, "evidence_id": "test:frequency",
                    "source": "synthetic", "source_version": "1", "population": "TEST",
                    "retrieved_at": "2026-09-14T00:00:00Z", "AC": ac, "AN": an,
                    "AF": ac / an, "quality_status": "PASS", "callable": True,
                    "homozygote_count": None, "hemizygote_count": None, **extra}

        value = self.evaluate(bs2, self.services([observation(hemizygote_count=1)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["zygosity_field"], "hemizygote_count")

    def test_bs2_x_linked_off_chrx_is_unknown(self):
        """A hemizygote count on an autosome would not be biologically meaningful."""
        self.input["inheritance"] = "x_linked_recessive"
        value = self.evaluate(bs2, self.services([self.observation(hemizygote_count=5)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)

    def test_bs2_no_genotype_count_data_is_unknown_not_not_met(self):
        """Every observation lacking the field is a search gap, not a negative result."""
        self.input["inheritance"] = "autosomal_recessive"
        value = self.evaluate(bs2, self.services([self.observation()]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["homozygote_count"])

    def test_bs2_incomplete_search_stays_unknown(self):
        self.input["inheritance"] = "autosomal_recessive"
        value = self.evaluate(
            bs2, self.services([self.observation(homozygote_count=0)], []))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("complete_population_evidence", value.missing_inputs)

    def test_bs2_reports_an_unconfigured_default_instead_of_inventing_one(self):
        del self.config["BS2"]["max_count"]
        self.input["inheritance"] = "autosomal_recessive"
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=1)]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["BS2.max_count"])

    def test_bs2_gene_specific_override_replaces_the_configured_default(self):
        """automated_core.context's gene_frequency_thresholds feeds
        input_data["bs2_threshold_override"], which criteria/bs2.py's threshold_policy()
        checks before config["BS2"] - a real curated, reviewed number is treated as
        APPROVED, not a draft prediction."""
        self.input["inheritance"] = "autosomal_recessive"
        self.input["bs2_threshold_override"] = {
            "max_count": 1, "source": "ClinGen test VCEP", "source_version": "1",
            "reviewed_at": "2026-09-25"}
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=2)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["threshold_scope"], "gene_specific")
        self.assertEqual(value.provenance["policy_status"], "APPROVED")
        self.assertNotIn("prediction pending curator sign-off", value.summary)
        self.assertEqual(value.review_points, [])

    def test_bs2_evidence_line_round_trips_through_va_spec(self):
        self.input["inheritance"] = "autosomal_recessive"
        value = self.evaluate(bs2, self.services([self.observation(homozygote_count=1)]))
        line = to_evidence_line(value)
        validate_1_0_1(line, "BS2")
        self.assertEqual(line["directionOfEvidenceProvided"], "disputes")


if __name__ == "__main__":
    unittest.main()

