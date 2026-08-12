# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
"""Wire-contract behavioral checks. Executed: see TEST_SPEC.md for the run record."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import struct
import unittest

# `fixtures` puts `src/` on sys.path, so this module runs standalone as well as under
# `python3 -m unittest discover -s tests`.
import fixtures

from spark_atomistic import (capabilities, run_atomistic_reference_unvalidated_json,
                             validate_atomistic_model_json)
from spark_atomistic.api import process_atomistic_json
from spark_atomistic.basin import DISABLED_CHECKPOINT_RECORD
from spark_atomistic.calculator import ProcessCalculator
from spark_atomistic.canonical import JSON_INTEGER_LIMIT, canonical_text, parse_json
from spark_atomistic.catalog import match_states, same_fixed_contract
from spark_atomistic.errors import DomainFailure, SEVERITY
from spark_atomistic.geometry import closest_periodic_displacement
from spark_atomistic.kinetics import RateTable
from spark_atomistic.model import (BASE_SPEC_SHA256, ERRATA_1_SHA256, ERRATA_2_SHA256,
                                   EXPECTED_SCHEMA_SHA256, SCHEMA_DESCRIPTOR,
                                   SCHEMA_REVISION, SCHEMA_SHA256, state_from_relaxation)
from spark_atomistic.resources import ResourceLedger
from spark_atomistic.rng import (ALGORITHM, PhiloxStream, derive_saddle_stream,
                                 derive_trajectory_stream, philox4x32_10, words_to_uniform)


ROOT = fixtures.ROOT
# E2-PAR-001 requires both backends to consume the SAME canonical corpus. It is kept as a
# single copy (currently under the Rust tree) and referenced from here rather than
# duplicated, because two copies of a golden file drift independently.
CORPUS = fixtures.CORPUS


class ContractSmokeTests(unittest.TestCase):
    def test_capabilities_are_conservative_and_immutable(self) -> None:
        # E2-API-007 freezes the capability value exactly. Keys such as `backend`,
        # `implementation_status`, `schema_revision`, and `schema_sha256` were part of the
        # pre-Errata-2 surface and MUST NOT reappear; `conformance` is `unvalidated`.
        value = capabilities()
        self.assertEqual(value["api"], "spark-atomistic-json/1")
        self.assertEqual(value["ir"], "spark-atomistic-model/1")
        self.assertEqual(value["basin_acceleration"], "disabled")
        self.assertEqual(value["conformance"], "unvalidated")
        self.assertFalse(value["validated"])
        self.assertFalse(value["production"])
        self.assertFalse(value["release"])
        self.assertEqual(list(value["operations"]), ["capabilities", "validate", "run"])
        self.assertEqual(dict(value["features"]), {
            "common_prefactor": True, "deterministic_checkpoint": True,
            "fixed_composition_off_lattice": True, "harmonic_tst": False,
            "local_environment_generic_reuse": False, "serial_kmc": True,
            "variable_composition": False,
        })
        self.assertEqual(set(value), {"api", "basin_acceleration", "conformance", "features",
                                      "ir", "operations", "production", "release", "validated"})
        self.assertEqual(len(SCHEMA_SHA256), 71)  # "sha256:" + 64 hexadecimal digits
        self.assertEqual(SCHEMA_REVISION, 2)
        with self.assertRaises(TypeError):
            value["api"] = "other"

    def test_spec_provenance_constants_are_machine_readable(self) -> None:
        # The three specification SHAs must be readable from code, not only from source
        # comments, so a provenance audit can be executed rather than eyeballed.
        self.assertEqual(SCHEMA_DESCRIPTOR["ir"], "spark-atomistic-model/1")
        self.assertEqual(SCHEMA_DESCRIPTOR["revision"], 2)
        for value in (BASE_SPEC_SHA256, ERRATA_1_SHA256, ERRATA_2_SHA256):
            self.assertEqual(len(value), 64)
            self.assertEqual(value, value.lower())
        # SCHEMA_DESCRIPTOR is hashed; adding Errata 2 to it would move SCHEMA_SHA256 and
        # break the cross-backend digest equality asserted by the parity corpus.
        self.assertNotIn("errata_2_sha256", SCHEMA_DESCRIPTOR)
        self.assertEqual(SCHEMA_SHA256, EXPECTED_SCHEMA_SHA256)

    def test_capability_response_matches_shared_parity_golden_byte_for_byte(self) -> None:
        # E2-PAR-001/E2-PAR-003: both backends consume the same canonical request corpus
        # and emit byte-identical responses. The golden files are the shared corpus.
        request = (CORPUS / "e2_capabilities.request.json").read_text().strip()
        golden = (CORPUS / "e2_capabilities.response.json").read_text().strip()
        self.assertEqual(process_atomistic_json(request), golden)

    def test_schema_fixture_and_explicit_run_gate(self) -> None:
        # E2-API-002/E2-API-003: the wire request is an envelope, not a bare model.
        model = json.loads((CORPUS / "e2_minimal_model.json").read_bytes())
        validated = json.loads(validate_atomistic_model_json(
            json.dumps({"model": model, "operation": "validate"})))
        self.assertEqual(validated["status"], "OK")
        self.assertEqual(validated["operation"], "validate")
        self.assertEqual(validated["value"]["ir"], "spark-atomistic-model/1")
        # E2-API-001: a bare model is not a request and is refused before any work.
        bare = json.loads(validate_atomistic_model_json(json.dumps(model)))
        self.assertEqual(bare["status"], "INVALID_INPUT")
        # E2-API-004: `run` without the boolean-true trust gate is refused.
        refused = json.loads(run_atomistic_reference_unvalidated_json(json.dumps(
            {"allow_unvalidated": False, "extension": {}, "model": model, "operation": "run"})))
        self.assertEqual(refused["status"], "INVALID_INPUT")
        self.assertEqual(refused["exit_code"], 64)

    def test_duplicate_keys_are_rejected(self) -> None:
        result = json.loads(validate_atomistic_model_json('{"schema":{},"schema":{}}'))
        self.assertEqual(result["status"], "INVALID_INPUT")
        # Pin the rejecting layer: duplicate keys must fail in the strict JSON reader
        # (E2-JSON-001), not merely as a side effect of failing the envelope shape check.
        self.assertEqual(result["context"]["requirement_id"], "E2-JSON-001")

    def test_exact_severity_strings(self) -> None:
        self.assertEqual(SEVERITY["DUPLICATE_EVENT"], "candidate reject")
        self.assertEqual(SEVERITY["CALCULATOR_FAILURE"], "transaction fail")
        self.assertEqual(SEVERITY["RATE_INVALID"], "fatal in strict mode")
        self.assertEqual(SEVERITY["NO_ENABLED_EVENT"],
                         "terminal-success if requested, else fatal")

    def test_unimplemented_retry_and_uncalibrated_alpha_are_rejected(self) -> None:
        # This gate is only meaningful against a baseline that otherwise validates. The
        # previous fixture (corpus/minimal_model.json) predates Errata 2 and is rejected
        # unconditionally, so every assertion below used to pass without exercising the
        # retry/alpha rules at all.
        def validate(model: dict) -> str:
            request = json.dumps({"model": model, "operation": "validate"})
            return json.loads(validate_atomistic_model_json(request))["status"]

        base = json.loads((CORPUS / "e2_minimal_model.json").read_bytes())
        self.assertEqual(validate(base), "OK")  # negative control

        retried = json.loads(json.dumps(base))
        retried["resources"]["retry_count"] = 1
        self.assertEqual(validate(retried), "INVALID_INPUT")

        alpha = json.loads(json.dumps(base))
        alpha["discovery"]["alpha"] = 0.1
        self.assertEqual(validate(alpha), "INVALID_INPUT")

    def test_errata_uniform_boundaries(self) -> None:
        low = words_to_uniform(0, 0)
        high = words_to_uniform(0xFFFFFFFF, 0xFFFFF000)
        self.assertEqual(struct.pack(">d", low).hex(), "3ca0000000000000")
        self.assertEqual(struct.pack(">d", high).hex(), "3fefffffffffffff")

    def test_triclinic_closest_vector_beats_component_rounding_case(self) -> None:
        cell = ((1.0, 0.0, 0.0), (0.9, 0.2, 0.0), (0.0, 0.0, 2.0))
        displacement, _ = closest_periodic_displacement((0.51, 0.11, 0.0), cell,
                                                        (True, True, False))
        brute = min(
            ((0.51 - i * cell[0][0] - j * cell[1][0]) ** 2
             + (0.11 - i * cell[0][1] - j * cell[1][1]) ** 2)
            for i in range(-4, 5) for j in range(-4, 5)
        )
        self.assertAlmostEqual(sum(item * item for item in displacement), brute, places=14)

    def test_state_identity_translation_image_and_same_species_permutation(self) -> None:
        request = {
            "atom_ids": ["a", "b"], "species": ["X", "X"],
            "positions": [(0.1, 0.2, 0.3), (0.7, 0.8, 0.9)],
            "cell": [(1.0, 0.0, 0.0), (0.2, 1.0, 0.0), (0.1, 0.3, 1.0)],
            "pbc": [True, True, True], "movable": [True, True],
            "constraints": {"kind": "fixed-mask"}, "charge": 0.0, "spin": 0.0,
            "calculator_model_digest": "fixture-model",
        }
        transformed = dict(request)
        transformed["atom_ids"] = ["b", "a"]
        transformed["species"] = ["X", "X"]
        transformed["movable"] = [True, True]
        transformed["positions"] = [(2.9, 2.3, 1.9), (1.3, 1.7, 1.3)]
        left = state_from_relaxation(request, -1.0, [(0.0, 0.0, 0.0)] * 2,
                                     1e-4, {"fixture": True})
        right = state_from_relaxation(transformed, -1.0, [(0.0, 0.0, 0.0)] * 2,
                                      1e-4, {"fixture": True})
        # E2-ID-002 states the certificate is mathematically invariant under whole-cell
        # translation, periodic-image choice, atom ID, and same-species permutation. In
        # binary64 the two anchor-minimum row arrays differ by one ULP (2.22e-16 here), so
        # the derived digests are NOT required to collide: E2-ID-005 makes candidate IDs
        # hints only and gives geometry/energy verification the deciding vote. Assert the
        # deciding path, and pin the hint's known non-strictness so it cannot silently
        # become a correctness assumption.
        self.assertTrue(same_fixed_contract(left, right))
        report = match_states(left, right, {
            "state_rms_tolerance": 1e-3, "state_max_tolerance": 1e-3,
            "state_energy_tolerance_per_atom": 1e-4,
        })
        self.assertTrue(report.equal)
        self.assertEqual(report.atom_mapping, (1, 0))
        self.assertLess(report.rms_displacement, 1e-12)
        self.assertNotEqual(left.candidate_identity, right.candidate_identity)

    def test_candidate_hash_never_blocks_tolerance_geometry_match(self) -> None:
        request = {
            "atom_ids": ["a", "b"], "species": ["X", "X"],
            "positions": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
            "cell": [(2.0, 0.0, 0.0), (0.2, 2.0, 0.0), (0.1, 0.3, 2.0)],
            "pbc": [True, True, True], "movable": [True, True],
            "constraints": {"kind": "fixed-mask"}, "charge": 0.0, "spin": 0.0,
            "calculator_model_digest": "fixture-model",
        }
        shifted = dict(request)
        shifted["positions"] = [(0.0, 0.0, 0.0), (0.5005, 0.5, 0.5)]
        left = state_from_relaxation(request, -1.0, [(0.0, 0.0, 0.0)] * 2,
                                     1e-4, {"fixture": True})
        right = state_from_relaxation(shifted, -1.0, [(0.0, 0.0, 0.0)] * 2,
                                      1e-4, {"fixture": True})
        self.assertTrue(same_fixed_contract(left, right))
        self.assertNotEqual(left.candidate_identity, right.candidate_identity)
        self.assertTrue(match_states(left, right, {
            "state_rms_tolerance": 1e-3, "state_max_tolerance": 5e-3,
            "state_energy_tolerance_per_atom": 1e-6,
        }).equal)

    def test_rate_table_snapshot_and_energy_status_split(self) -> None:
        # E2-KMC-002: the envelope has exactly `payload` and `payload_sha256`.
        snapshot = RateTable("s", ("e",), ("d",), (0.0,), (1.0,), 1.0, None).snapshot()
        self.assertEqual(set(snapshot), {"payload", "payload_sha256"})
        self.assertEqual(snapshot["payload"]["origin_state_id"], "s")
        self.assertEqual(snapshot["payload"]["schema"],
                         "spark-atomistic-rate-table-snapshot/1")
        self.assertEqual(set(snapshot["payload"]),
                         {"destination_state_ids", "event_ids", "log_rates",
                          "lost_rate_log_upper_bound", "origin_state_id", "rates",
                          "schema", "total_rate_per_s"})
        self.assertTrue(snapshot["payload_sha256"].startswith("sha256:"))
        self.assertEqual(len(snapshot["payload_sha256"]), 71)
        ledger = ResourceLedger(10, 10.0, 1 << 40, 1 << 20, 10, 0.0)
        calculator = ProcessCalculator(
            {"command": ["/unused"], "model_name": "m", "model_version": "1",
             "model_digest": "d", "deterministic": True}, 1.0, 100, 100, ledger)
        response = {"status": "OK", "energy": "bad", "forces": [[0.0, 0.0, 0.0]],
                    "units": {"energy": "eV", "forces": "eV/angstrom"},
                    "model_name": "m", "model_version": "1", "model_digest": "d",
                    "evaluation_id": "x", "deterministic": True, "request_digest": "q"}
        with self.assertRaises(DomainFailure) as malformed:
            calculator._validate_response(response, {"atom_ids": ["a"]}, "q", "test", "x")
        self.assertEqual(malformed.exception.outcome.status, "CALCULATOR_FAILURE")
        response["energy"] = float("nan")
        with self.assertRaises(DomainFailure) as nonfinite:
            calculator._validate_response(response, {"atom_ids": ["a"]}, "q", "test", "x")
        self.assertEqual(nonfinite.exception.outcome.status, "NONFINITE_RESULT")


class ParityCorpusTests(unittest.TestCase):
    """Mandatory `E2-PAR-002` cases that are pure functions of the wire contract."""

    def test_canonical_number_cases_match_the_shared_golden(self) -> None:
        # E2-CAN-004 / E2-PAR-002 item 4: RFC 8785 shortest round-tripping decimal,
        # negative zero as `0`, and the four values adjacent to 1e-6 and 1e21.
        golden = json.loads((CORPUS / "e2_canonical_numbers.json").read_bytes())
        self.assertEqual(len(golden["cases"]), 7)
        for case in golden["cases"]:
            with self.subTest(case=case["binary64"]):
                if case["binary64"] == "negative-zero":
                    value = -0.0
                else:
                    value = struct.unpack(">d", bytes.fromhex(case["binary64"][2:]))[0]
                self.assertEqual(canonical_text(value), case["canonical"])
                if case["binary64"] != "negative-zero":
                    self.assertEqual(float(case["canonical"]), value)
        self.assertEqual(golden["integer_max"], JSON_INTEGER_LIMIT)
        self.assertEqual(golden["integer_min"], -JSON_INTEGER_LIMIT)

    def test_portable_integer_boundary_and_out_of_domain(self) -> None:
        # E2-JSON-002: the boundary itself is legal; one past it is INVALID_INPUT and is
        # rejected before schema validation.
        self.assertEqual(canonical_text(JSON_INTEGER_LIMIT), "9007199254740991")
        self.assertEqual(canonical_text(-JSON_INTEGER_LIMIT), "-9007199254740991")
        accepted = json.loads(validate_atomistic_model_json(
            '{"model":{"schema":{"id":"spark-atomistic-model/1"}},"operation":"validate"}'))
        self.assertEqual(accepted["status"], "INVALID_INPUT")  # incomplete model, not integers
        rejected = json.loads(validate_atomistic_model_json(
            '{"model":{"n":9007199254740992},"operation":"validate"}'))
        self.assertEqual(rejected["status"], "INVALID_INPUT")
        self.assertEqual(rejected["context"]["requirement_id"], "E2-JSON-002")

    def test_nonfinite_tokens_and_lone_surrogates_are_refused(self) -> None:
        # E2-JSON-003 nonfinite; E2-JSON-001 lone UTF-16 surrogate and malformed UTF-8.
        nonfinite = json.loads(validate_atomistic_model_json(
            '{"model":{"x":NaN},"operation":"validate"}'))
        self.assertEqual(nonfinite["status"], "NONFINITE_RESULT")
        self.assertEqual(nonfinite["exit_code"], 65)
        surrogate = json.loads(validate_atomistic_model_json(
            '{"model":{"x":"\\ud800"},"operation":"validate"}'))
        self.assertEqual(surrogate["status"], "INVALID_INPUT")
        self.assertEqual(surrogate["context"]["requirement_id"], "E2-JSON-001")
        malformed = json.loads(validate_atomistic_model_json(b'{"operation":"validate","model":{"x":"\xff\xfe"}}'))
        self.assertEqual(malformed["status"], "INVALID_INPUT")
        self.assertEqual(malformed["exit_code"], 64)

    def test_duplicate_key_golden_file_is_refused_at_any_depth(self) -> None:
        # E2-JSON-001: duplicate keys at ANY nesting depth. The golden file nests one.
        raw = (CORPUS / "strict_duplicate.invalid.json").read_bytes()
        with self.assertRaises(DomainFailure) as failure:
            parse_json(raw)
        self.assertEqual(failure.exception.outcome.status, "INVALID_INPUT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-JSON-001")

    def test_metadata_is_free_form_and_excluded_from_the_config_digest(self) -> None:
        # E2-SCHEMA-001 / E2-CAN-007 / E2-PAR-002 item 2: metadata-only variants share
        # one config digest. Paired with a behavioural field change, which must NOT.
        base = json.loads((CORPUS / "e2_minimal_model.json").read_bytes())

        def config_digest(model):
            response = json.loads(validate_atomistic_model_json(
                json.dumps({"model": model, "operation": "validate"})))
            self.assertEqual(response["status"], "OK")
            return response["value"]["config_digest"]

        plain = config_digest(base)
        for metadata in ({}, {"note": "a"}, {"nested": {"deep": [1, 2.5, None, True]}}):
            with self.subTest(metadata=metadata):
                variant = json.loads(json.dumps(base))
                variant["metadata"] = metadata
                self.assertEqual(config_digest(variant), plain)
        behavioural = json.loads(json.dumps(base))
        behavioural["kinetics"]["temperature"] = 301.0
        self.assertNotEqual(config_digest(behavioural), plain)

    def test_unknown_keys_are_refused_in_every_normative_object(self) -> None:
        # E2-SCOPE-003: only root `metadata` and the run-request `extension` are free.
        base = json.loads((CORPUS / "e2_minimal_model.json").read_bytes())
        for path in ("kinetics", "resources", "output", "discovery", "saddle_search"):
            with self.subTest(path=path):
                variant = json.loads(json.dumps(base))
                variant[path]["unexpected"] = 1
                response = json.loads(validate_atomistic_model_json(
                    json.dumps({"model": variant, "operation": "validate"})))
                self.assertEqual(response["status"], "INVALID_INPUT")
                self.assertEqual(response["context"]["requirement_id"], "E2-SCOPE-003")

    def test_basin_enabled_is_schema_valid_but_never_advertised(self) -> None:
        # E2-SCHEMA-011 / E2-BASIN-001: `true` validates OK; capability stays `disabled`.
        base = json.loads((CORPUS / "e2_minimal_model.json").read_bytes())
        base["basin"]["enabled"] = True
        response = json.loads(validate_atomistic_model_json(
            json.dumps({"model": base, "operation": "validate"})))
        self.assertEqual(response["status"], "OK")
        self.assertEqual(capabilities()["basin_acceleration"], "disabled")
        self.assertEqual(DISABLED_CHECKPOINT_RECORD, {"enabled": False, "reason": "v1-disabled"})


class PhiloxCorpusTests(unittest.TestCase):
    """`E2-RNG-001..005` derivation, counter arithmetic, and state restore."""

    def setUp(self) -> None:
        self.golden = json.loads((CORPUS / "philox_errata1.json").read_bytes())

    def test_zero_key_zero_counter_block_matches_the_shared_golden(self) -> None:
        self.assertEqual(list(philox4x32_10((0, 0, 0, 0), (0, 0))),
                         self.golden["zero_counter_zero_key_words"])
        self.assertEqual(ALGORITHM, self.golden["algorithm"])

    def test_boundary_records_match_errata_1_and_their_hashes(self) -> None:
        # E2-RNG-006 / E2-RNG-007: both golden lines and their SHA-256 are normative.
        for case in self.golden["boundary_cases"]:
            with self.subTest(q=case["q"]):
                value = words_to_uniform(case["a"], case["b"])
                self.assertEqual("0x" + struct.pack(">d", value).hex(), case["bits"])
                self.assertEqual((case["a"] << 20) | (case["b"] >> 12), case["q"])
                self.assertTrue(0.0 < value < 1.0)
        low = {"a": 0, "b": 0, "q": 0, "raw_binary64_bits": "0x3ca0000000000000",
               "uniform_hex": "0x1.0000000000000p-53"}
        high = {"a": 4294967295, "b": 4294963200, "q": 4503599627370495,
                "raw_binary64_bits": "0x3fefffffffffffff",
                "uniform_hex": "0x1.fffffffffffffp-1"}
        self.assertEqual(
            hashlib.sha256(canonical_text(low).encode()).hexdigest(),
            "6ce1fb5214530ba6b04e4bf75aaeba5d02acf6694cd462004faf7640a665fc03")
        self.assertEqual(
            hashlib.sha256(canonical_text(high).encode()).hexdigest(),
            "a15157604d319e5525e3b83eba02259088e317ea2b0b0e1a3bb28060e093cf43")

    def test_counter_carry_propagates_towards_c3(self) -> None:
        # E2-RNG-003: c0 is least significant; the block increment is a 128-bit add.
        stream = PhiloxStream((0, 0), (1 << 32) - 1, (1 << 32) - 1)
        stream.uniform()
        stream.uniform()
        self.assertEqual(stream.checkpoint()["next_counter"], [0, 1, 0, 0])
        self.assertEqual(stream.consumed_blocks, 1)

    def test_block_pairing_and_buffer_state_machine(self) -> None:
        # E2-RNG-002: first pair keeps the block with next_pair=1; second releases it.
        stream = derive_trajectory_stream(0)
        self.assertIsNone(stream.checkpoint()["buffered_block"])
        stream.uniform()
        first = stream.checkpoint()
        self.assertIsNotNone(first["buffered_block"])
        self.assertEqual(first["next_pair"], 1)
        self.assertEqual(first["consumed_blocks"], 1)
        stream.uniform()
        second = stream.checkpoint()
        self.assertIsNone(second["buffered_block"])
        self.assertEqual(second["next_pair"], 0)
        self.assertEqual(second["consumed_uniforms"], 2)
        # Round trip is exact, and a mid-block state restores its retained block.
        self.assertEqual(PhiloxStream.restore(first).checkpoint(), first)
        self.assertEqual(PhiloxStream.restore(second).checkpoint(), second)

    def test_substream_derivation_is_scheduling_independent(self) -> None:
        # DET-003 / E2-RNG-004: a substream depends only on seed, state ID, class, index.
        left = derive_saddle_stream(7, "state:sha256:aa", "global", 3)
        right = derive_saddle_stream(7, "state:sha256:aa", "global", 3)
        self.assertEqual(left.checkpoint(), right.checkpoint())
        for changed in (derive_saddle_stream(8, "state:sha256:aa", "global", 3),
                        derive_saddle_stream(7, "state:sha256:ab", "global", 3),
                        derive_saddle_stream(7, "state:sha256:aa", "local", 3),
                        derive_saddle_stream(7, "state:sha256:aa", "global", 4),
                        derive_saddle_stream(7, "state:sha256:aa", "class-selection", 3)):
            self.assertNotEqual(left.checkpoint()["key"] + left.checkpoint()["initial_counter"],
                                changed.checkpoint()["key"] + changed.checkpoint()["initial_counter"])
        # The length prefixes stop a class/state boundary from being ambiguous.
        self.assertNotEqual(derive_saddle_stream(0, "ab", "c", 0).checkpoint()["key"],
                            derive_saddle_stream(0, "a", "bc", 0).checkpoint()["key"])
        # The trajectory stream is a different derivation material entirely.
        self.assertNotEqual(derive_trajectory_stream(7).checkpoint()["key"],
                            left.checkpoint()["key"])



if __name__ == "__main__":
    unittest.main()
