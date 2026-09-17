"""API result persistence: what gets written, and what must keep working when it cannot be.

The evaluation itself is stubbed out here. These tests are about the artifact the API
leaves behind - whether one is written at all, what it contains, and that a failure to
write never turns a completed judgment into an error response.
"""

import json
import os
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from api import main as api_main  # noqa: E402
from api.result_store import OUTPUT_DIR_ENV, filename, output_dir, save  # noqa: E402


VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "14\t23420189\tcase3-noise6\tC\tT\t80\tPASS\t"
    "ASSEMBLY=GRCh38;GENE=MYH7;TRANSCRIPT=NM_000257.4;HGVSC=c.3382G>A;HGVSP=p.(Ala1128Thr)\n"
)
LINES = {"BA1": {"id": "evline:MYH7_c_3382G_A_BA1", "type": "EvidenceLine"}}


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv(OUTPUT_DIR_ENV, str(tmp_path))
    return TestClient(api_main.app), tmp_path


def post(client, criteria=("BA1",)):
    with patch.object(api_main, "run_selected_criteria", new=AsyncMock(return_value=LINES)):
        return client.post(
            "/v1/get_evidence_line_by_target_criteria",
            json={"vcf": VCF, "clinical_note": "", "criteria": list(criteria)},
        )


def written(directory):
    return sorted(directory.glob("*.json"))


def test_a_synchronous_result_is_written(client):
    api, directory = client
    response = post(api)
    assert response.status_code == 200
    files = written(directory)
    assert len(files) == 1
    document = json.loads(files[0].read_text(encoding="utf-8"))
    assert document["status"] == "succeeded"
    assert document["result"] == LINES
    assert document["endpoint"] == "get_evidence_line_by_target_criteria"


def test_the_artifact_identifies_the_request_that_produced_it(client):
    """The payload alone does not say which variant or criteria were asked for."""
    api, directory = client
    post(api, criteria=("BA1", "PM2"))
    document = json.loads(written(directory)[0].read_text(encoding="utf-8"))
    assert document["variant"]["gene"] == "MYH7"
    assert document["variant"]["hgvsc"] == "c.3382G>A"
    assert document["criteria"] == ["PM2", "BA1"]  # canonical ACMG order, as evaluated
    assert document["finished_at"]


def test_nothing_is_written_when_persistence_is_not_configured(tmp_path, monkeypatch):
    monkeypatch.delenv(OUTPUT_DIR_ENV, raising=False)
    assert output_dir() is None
    response = post(TestClient(api_main.app))
    assert response.status_code == 200
    assert written(tmp_path) == []


def test_a_blank_setting_counts_as_unconfigured(monkeypatch):
    monkeypatch.setenv(OUTPUT_DIR_ENV, "   ")
    assert output_dir() is None


def test_the_directory_is_created_on_demand(tmp_path, monkeypatch):
    target = tmp_path / "nested" / "results"
    monkeypatch.setenv(OUTPUT_DIR_ENV, str(target))
    post(TestClient(api_main.app))
    assert len(written(target)) == 1


def test_an_unwritable_directory_does_not_fail_the_request(tmp_path, monkeypatch, caplog):
    """A full disk must not turn a completed evaluation into an error response."""
    blocked = tmp_path / "blocked"
    blocked.write_text("not a directory", encoding="utf-8")
    monkeypatch.setenv(OUTPUT_DIR_ENV, str(blocked))
    response = post(TestClient(api_main.app))
    assert response.status_code == 200
    assert response.json()["evidence_lines"] == LINES
    assert any("Could not write API result" in record.message for record in caplog.records)


def test_two_results_for_one_variant_do_not_overwrite_each_other(client):
    api, directory = client
    post(api)
    post(api)
    assert len(written(directory)) == 2


def test_a_filename_survives_characters_hgvs_uses_and_filesystems_reject(tmp_path, monkeypatch):
    monkeypatch.setenv(OUTPUT_DIR_ENV, str(tmp_path))
    name = filename(
        endpoint="get_evidence_line_by_target_criteria",
        variant={"gene": "MYH7", "hgvsc": "c.3382G>A"},
        identifier="0123456789abcdef",
        finished_at="2026-09-17T02:40:13.990000+00:00",
    )
    assert not set(name) & set('<>:"/\\|?*')
    assert "MYH7" in name and "3382G_A" in name
    (tmp_path / name).write_text("{}", encoding="utf-8")  # the filesystem accepts it


def test_a_variant_without_gene_or_hgvsc_still_gets_a_name():
    name = filename(endpoint="classify_criteria", variant={}, identifier="abcd1234",
                    finished_at="2026-09-17T02:40:13+00:00")
    assert "nogene" in name and "novariant" in name


def test_an_unserializable_payload_is_reported_not_raised(tmp_path, monkeypatch, caplog):
    monkeypatch.setenv(OUTPUT_DIR_ENV, str(tmp_path))
    assert save(endpoint="classify_criteria", variant={}, criteria=None,
                payload={"bad": object()}, status="succeeded", identifier="x") is None
    assert written(tmp_path) == []
    assert any("Could not write API result" in record.message for record in caplog.records)


def test_a_failed_job_is_recorded_too(tmp_path, monkeypatch):
    """A failure that a restart would erase is worth keeping as much as a success."""
    monkeypatch.setenv(OUTPUT_DIR_ENV, str(tmp_path))
    path = save(endpoint="classify_criteria", variant={"gene": "MYH7", "hgvsc": "c.1A>G"},
                criteria=None, payload={"type": "RuntimeError", "message": "vLLM unreachable"},
                status="failed", identifier="job-1")
    document = json.loads(Path(path).read_text(encoding="utf-8"))
    assert document["status"] == "failed"
    assert document["result"]["type"] == "RuntimeError"
