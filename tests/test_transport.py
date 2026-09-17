import gzip
import io
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from acmg_pipeline.providers.http import CachedHttpClient, FetchError
from acmg_pipeline.providers.vep import VepProvider
from acmg_pipeline.automated_core.models import Variant


class MemoryCache:
    """No filesystem temp-directory dependency for transport contract tests."""
    def __init__(self):
        self.data = {}

    def exists(self, path):
        return str(path) in self.data

    def read(self, path, **kwargs):
        return self.data[str(path)]

    def write(self, path, value, **kwargs):
        self.data[str(path)] = value


class TransportTests(unittest.TestCase):
    def test_cache_replay_and_integrity(self):
        cache = MemoryCache()
        calls = []
        def opener(request, timeout):
            calls.append(request.full_url)
            return io.BytesIO(b'{"value": 2}')
        with patch.object(Path, "exists", lambda p: cache.exists(p)), \
             patch.object(Path, "read_text", lambda p, **kw: cache.read(p, **kw)), \
             patch.object(Path, "write_text", lambda p, v, **kw: cache.write(p, v, **kw)), \
             patch.object(Path, "mkdir"):
            client = CachedHttpClient("cache", opener=opener, sleeper=lambda _: None)
            first = client.fetch("https://example.org/data", dataset_version="1")
            client.offline = True
            self.assertEqual(first, client.fetch("https://example.org/data", dataset_version="1"))
            self.assertEqual(len(calls), 1)
            with self.assertRaisesRegex(FetchError, "OFFLINE_CACHE_MISS"):
                client.fetch("https://example.org/data", dataset_version="2")
            path = next(iter(cache.data))
            document = json.loads(cache.data[path])
            document["body"]["value"] = 3
            cache.data[path] = json.dumps(document)
            with self.assertRaisesRegex(FetchError, "CACHE_INTEGRITY_FAILURE"):
                client.fetch("https://example.org/data", dataset_version="1")

    def test_gzip_text_is_decompressed_once_and_cached_as_text(self):
        """NCBI publishes the MANE summary only as .gz, and decoding a gzip stream as UTF-8
        first destroys it. The cache holds the decompressed text, so a replay never has to
        decompress anything."""
        cache = MemoryCache()
        calls = []

        def opener(request, timeout):
            calls.append(request.full_url)
            return io.BytesIO(gzip.compress(
                b"symbol\tRefSeq_nuc\nMYBPC3\tNM_000256.3\n"))

        with patch.object(Path, "exists", lambda p: cache.exists(p)), \
             patch.object(Path, "read_text", lambda p, **kw: cache.read(p, **kw)), \
             patch.object(Path, "write_text", lambda p, v, **kw: cache.write(p, v, **kw)), \
             patch.object(Path, "mkdir"):
            client = CachedHttpClient("cache", opener=opener, sleeper=lambda _: None)
            first = client.fetch("https://example.org/summary.txt.gz", response_format="text-gz")
            self.assertIn("MYBPC3", first["body"])
            client.offline = True
            replay = client.fetch("https://example.org/summary.txt.gz", response_format="text-gz")
            self.assertEqual(first["body"], replay["body"])
            self.assertEqual(len(calls), 1)

    def test_an_unknown_response_format_is_refused(self):
        client = CachedHttpClient("absent-cache", opener=None, sleeper=lambda _: None)
        with self.assertRaises(ValueError):
            client.fetch("https://example.org/x", response_format="xml")

    def test_retry_is_bounded(self):
        calls = []
        def opener(request, timeout):
            calls.append(request)
            raise TimeoutError("test")
        client = CachedHttpClient("absent-cache", opener=opener, sleeper=lambda _: None, attempts=2)
        with self.assertRaises(FetchError):
            client.fetch("https://example.org/retry")
        self.assertEqual(len(calls), 2)

    def test_vep_contract(self):
        class Client:
            def fetch(self, url, **kwargs):
                return {"retrieved_at": "2026-09-15", "body": [{"transcript_consequences": [
                    {"transcript_id": "ENST_TEST.1", "gene_symbol": "TEST", "consequence_terms": ["missense_variant"],
                     "amino_acids": "R/W", "protein_start": 10, "protein_end": 10}]}]}
        result = VepProvider(Client(), "test").annotate(Variant("GRCh38", "1", 2, "C", "T"))
        self.assertEqual(result[0]["ref_aa"], "R")
        self.assertEqual(result[0]["alt_aa"], "W")
        self.assertEqual(result[0]["source_version"], "test")
