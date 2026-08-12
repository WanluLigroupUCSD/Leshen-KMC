# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Python half of the `E2-PAR-005` cross-language experiment.

`E2-PAR-005`: "Execution conformance additionally requires every fixture to pass in every
implementation." A single-backend suite cannot establish that, so these tests do three things
this package could not do before: they run the Python backend against the SHARED corpus that
drives the Rust suite, they prove the Python backend re-reads every artifact it emits, and they
compare against the Rust artifacts byte-for-byte when those artifacts are present.

The comparison is driven by `tests/xlang_harness.py`; the Rust artifacts come from
`cargo test --test xlang_emit` with `SPARK_XLANG_OUT` set to the same directory. When they are
absent the comparison test SKIPS with the exact command that produces them, so an unavailable
comparison is never reported as a passing one.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest

import fixtures  # noqa: F401  (puts src/ on sys.path)
import xlang_harness

from spark_atomistic.api import process_atomistic_json
from spark_atomistic.canonical import canonical_bytes, parse_json


CORPUS = xlang_harness.XLANG


def response(raw):
    return json.loads(process_atomistic_json(raw))


def metadata_request(raw_metadata_text):
    model = (CORPUS.parent / "e2_minimal_model.json").read_bytes()
    body = canonical_bytes(json.loads(model)).decode()[:-1]
    return ('{"model":' + body + ',"metadata":' + raw_metadata_text
            + '},"operation":"validate"}').encode()


class PortableIntegerDomain(unittest.TestCase):
    """`E2-JSON-002` in the Python backend, checked against the exact D-E2-01 hole.

    D-E2-01 was a Rust defect: the domain was enforced only in the strict parser's i64/u64
    visitors, so every integer literal wider than u64 reached the binary64 visitor and was
    accepted, silently rounded. Python decodes integer literals as arbitrary-precision `int`
    and applies `validate_portable_value` to the parsed tree, so the same inputs must be
    rejected here. Each rejection is paired with an accepted baseline that differs only in the
    property under test.
    """

    def test_domain_endpoints_are_accepted(self):
        for token in ("9007199254740991", "-9007199254740991", "0"):
            with self.subTest(token=token):
                self.assertEqual(response(metadata_request('{"x":%s}' % token))["status"], "OK")

    def test_every_out_of_domain_integer_literal_is_rejected(self):
        # 2^53 and 2^64-1 were already rejected by Rust; 2^64, a digit-form 1e30 and a
        # 401-digit literal are exactly the values that escaped it.
        for token in ("9007199254740992", "-9007199254740992", "18446744073709551615",
                      "18446744073709551616", "-18446744073709551616", "1" + "0" * 30,
                      "999999999999999900000", "1" + "0" * 400):
            with self.subTest(token=token[:24]):
                value = response(metadata_request('{"x":%s}' % token))
                self.assertEqual(value["status"], "INVALID_INPUT")
                self.assertEqual(value["exit_code"], 64)
                self.assertEqual(value["message"], "input invalid")
                self.assertEqual(value["context"]["requirement_id"], "E2-JSON-002")

    def test_the_same_magnitudes_with_an_exponent_are_binary64_not_integers(self):
        # ACCEPTED BASELINES. E2-JSON-002 speaks of a "syntactically valid integer"; a token
        # carrying a fraction or an exponent is an E2-JSON-003 binary64 and stays admissible.
        for token in ("1e30", "-1e30", "1.8446744073709552e19", "9.999999999999999e20"):
            with self.subTest(token=token):
                self.assertEqual(response(metadata_request('{"x":%s}' % token))["status"], "OK")

    def test_digits_inside_a_string_are_not_a_number_token(self):
        self.assertEqual(response(metadata_request('{"x":"18446744073709551616"}'))["status"], "OK")

    def test_e3_can_001_keeps_large_binary64_output_in_the_real_domain(self):
        """E3-CAN-001 closes the encoder/parser domain under round-trip."""
        emitted = canonical_bytes(9.999999999999999e20)
        self.assertEqual(emitted, b"9.999999999999999e+20")
        self.assertEqual(response(metadata_request('{"x":%s}' % emitted.decode()))["status"], "OK")
        self.assertEqual(parse_json(emitted), 9.999999999999999e20)
        self.assertEqual(canonical_bytes(float(2**53)), b"9.007199254740992e+15")
        self.assertEqual(canonical_bytes(1e21), b"1e+21")
        self.assertEqual(response(metadata_request('{"x":1e+21}'))["status"], "OK")


class SharedCorpusExecution(unittest.TestCase):
    """`E2-PAR-001`: the Python backend consumes the same canonical corpus as the Rust backend."""

    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.mkdtemp(prefix="spark-xlang-py-")
        cls.ids = xlang_harness.emit(cls.directory)
        cls.out = Path(cls.directory) / "python"

    def test_every_manifest_case_emits_exactly_one_artifact(self):
        manifest = xlang_harness.load_manifest()
        probes = json.loads((CORPUS / "probes.json").read_bytes())
        expected = (len(manifest["requests"]) + len(probes["cases"])
                    + len(manifest["checkpoints"]))
        self.assertEqual(len(self.ids), expected)
        self.assertEqual(len(set(self.ids)), len(self.ids))
        for case_id in self.ids:
            self.assertTrue((self.out / (case_id + ".out")).exists(), case_id)

    def test_the_backend_can_reread_every_artifact_it_emits(self):
        """`E2-CAN-001` + `E2-CKPT-007`(1): canonical bytes must re-parse and re-canonicalise.

        Rust historically failed four cases under D-E2-03. `float_roundtrip` repaired that
        defect; both backends now require zero self-read failures.
        """
        broken = []
        for case_id in self.ids:
            raw = (self.out / (case_id + ".out")).read_bytes()
            if canonical_bytes(parse_json(raw)) != raw:
                broken.append(case_id)
        self.assertEqual(broken, [])

    def test_the_shared_checkpoint_with_three_reciprocal_pairs_restores(self):
        """D-E2-02 counter-example: this catalog IS representable and IS restorable.

        `checkpoints/checkpoint-clean.json` holds six directed records (three reciprocal
        pairs), four committed states and one committed KMC step. Errata 3 makes it
        representable and restorable in both backends. The paired control is the same model
        with an empty catalog.
        """
        clean = json.loads((self.out / "checkpoint-clean__restore.out").read_bytes())
        self.assertEqual(clean["outcome"], "restored")
        self.assertEqual(len(clean["catalog_event_ids"]), 6)
        self.assertEqual(len(clean["catalog_state_ids"]), 4)
        self.assertEqual(clean["step_index"], 1)
        control = json.loads((self.out / "checkpoint-clean__zero_step.out").read_bytes())
        self.assertEqual(control["outcome"], "restored")
        self.assertEqual(control["catalog_event_ids"], [])


class RustComparison(unittest.TestCase):
    """`E2-PAR-003`: the two backends' canonical bytes must be identical."""

    def test_recorded_divergences_are_exactly_the_documented_set(self):
        root = os.environ.get("SPARK_XLANG_OUT")
        if not root or not (Path(root) / "rust" / "_index.json").exists():
            raise unittest.SkipTest(
                "Rust artifacts absent. Produce them with:\n"
                "  cd ../spark-atomistic-rs && SPARK_XLANG_OUT=<dir> "
                "CARGO_TARGET_DIR=<tmp> cargo test --test xlang_emit\n"
                "then rerun with SPARK_XLANG_OUT=<dir>. This test SKIPS rather than passes so "
                "that an unavailable comparison is never mistaken for a clean one.")
        xlang_harness.emit(root)
        xlang_harness.compare(root)
        report = json.loads((Path(root) / "parity_report.json").read_bytes())
        divergent = sorted(r["case"] for r in report["cases"] if r["verdict"] == "DIVERGENT")
        self.assertEqual(report["summary"]["NOT-APPLICABLE"], 0)
        self.assertEqual(divergent, sorted(DOCUMENTED_DIVERGENCES),
                         "cross-language divergence set changed; update TEST_SPEC.md before "
                         "editing this list")


# Finding F6 preserves the historical pre-Errata-3 divergence record. D-127 adopted
# Errata 3; the current normative byte-identity gate permits no divergence.
DOCUMENTED_DIVERGENCES = []


if __name__ == "__main__":
    unittest.main()
