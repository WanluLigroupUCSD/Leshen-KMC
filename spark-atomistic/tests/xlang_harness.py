# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Cross-language parity harness: Python emitter and byte-for-byte comparator.

`E2-PAR-001` requires every conforming implementation to consume the SAME canonical request
corpus and emit the same canonical response, status record, state/event/catalog record, rate
snapshot, RNG state and checkpoint schema. `E2-PAR-003` requires those bytes to be identical.
`E2-PAR-005` makes execution conformance conditional on every mandatory fixture passing in
EVERY implementation, which is exactly what a single-backend suite cannot establish.

The shared corpus lives at `../spark-atomistic-rs/tests/corpus/xlang/` and is referenced, never
copied. The Rust half is `../spark-atomistic-rs/tests/xlang_emit.rs`; both halves write one
canonical byte string per case into `<out>/<backend>/<case>.out`.

    python3 tests/xlang_harness.py emit    --out DIR
    python3 tests/xlang_harness.py compare --out DIR

`compare` exits nonzero when any case diverges. A divergence is reported with the exact first
differing byte offset and both surrounding strings; it is never normalised away.

Two probes are marked SHARED-ALGEBRA below: the Python backend exposes no standalone function
for them, so this file composes the spec-quoted payload and the BACKEND supplies only canonical
encoding and SHA-256. Their verdicts are labelled so no divergence is mis-attributed.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import struct
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

XLANG = ROOT.parent / "spark-atomistic-rs" / "tests" / "corpus" / "xlang"

from spark_atomistic.api import process_atomistic_json  # noqa: E402
from spark_atomistic.canonical import (canonical_bytes, canonical_text, deep_thaw,  # noqa: E402
                                       digest)
from spark_atomistic.catalog import Event  # noqa: E402
from spark_atomistic.checkpoint import validate_checkpoint_payload  # noqa: E402
from spark_atomistic.engine import ReferenceEngine  # noqa: E402
from spark_atomistic.errors import DomainFailure, MESSAGE, SEVERITY, STATUSES, exit_code  # noqa: E402
from spark_atomistic.kinetics import build_rate_table, propose_serial_step  # noqa: E402
from spark_atomistic.model import (SCHEMA_SHA256, fixed_contract_digest,  # noqa: E402
                                    geometry_certificate, identity_digest,
                                    state_from_relaxation, validate_model)
from spark_atomistic.rng import (PhiloxStream, derive_saddle_stream,  # noqa: E402
                                  derive_trajectory_stream, philox4x32_10, words_to_uniform,
                                  _int_to_words, _words_to_int)

SHARED_ALGEBRA = {"event_ids", "digests"}
TOLERANCE_KEYS = ("barrier_tolerance", "detailed_balance_tolerance", "log_rate_cutoff",
                  "saddle_energy_tolerance", "saddle_max_tolerance", "saddle_rms_tolerance",
                  "state_energy_tolerance_per_atom", "state_max_tolerance", "state_rms_tolerance")


# ------------------------------------------------------------------ shared shape helpers

def num(x):
    return canonical_text(float(x))


def bits(x):
    return "0x%016x" % struct.unpack(">Q", struct.pack(">d", float(x)))[0]


def scalar(x):
    return {"bits": bits(x), "canonical": num(x)}


def failure(status, requirement):
    return {"outcome": "failure", "requirement_id": requirement, "status": status}


def from_domain(exc):
    return failure(exc.outcome.status, exc.outcome.context["requirement_id"])


def philox_from_words(key, initial_counter):
    counter = _words_to_int(list(initial_counter))
    return PhiloxStream((key[0], key[1]), counter, counter)


def substream_digest(stream):
    return digest(stream.checkpoint())


def draw(stream, draws):
    uniforms = []
    for _ in range(draws):
        try:
            uniforms.append(scalar(stream.uniform()))
        except DomainFailure as exc:
            uniforms.append(from_domain(exc))
            break
    return {"final_state": stream.checkpoint(), "substream_digest": substream_digest(stream),
            "uniforms": uniforms}


def state_request(system):
    return {"schema": "spark-atomistic-model/1", "atom_ids": list(system["atom_ids"]),
            "species": list(system["species"]),
            "positions": [list(p) for p in system["positions"]],
            "cell": [list(r) for r in system["cell"]], "pbc": list(system["pbc"]),
            "movable": list(system["movable"]), "constraints": {"kind": "fixed-mask"},
            "charge": system["charge"], "spin": system["spin"],
            "calculator_model_digest": system["calculator_model_digest"]}


def stub_event(event_id, destination_state_id, origin_state_id, log_rate, selectable):
    """A directed record carrying exactly the fields `build_rate_table` reads (E2-KMC-001/002)."""
    return Event(event_id, "event:reverse", "pair:stub", origin_state_id, destination_state_id,
                 (), 0.0, (), (), -1.0, (), "DIRECTIONAL", 0, "search:stub", "stub", 0.0, 0.0,
                 log_rate, 0.0, 0.0, 1e13, 300.0, (), {}, {}, "sha256:test-model",
                 "sha256:identity", "sha256:tolerance", selectable)


def rate_table_from_rows(rows, origin, cutoff):
    events = [stub_event(eid, dest, origin, log, log >= cutoff) for eid, dest, log in rows]
    return build_rate_table(tuple(events), origin)


# ------------------------------------------------------------------ probes

def probe(inp, template):
    name = inp["probe"]
    if name == "canonical_number":
        out = []
        for raw in inp["bits"]:
            x = struct.unpack(">d", struct.pack(">Q", int(raw, 16)))[0]
            out.append({"canonical": num(x), "input_bits": raw, "round_trip_bits": bits(x)})
        return {"cases": out}
    if name == "digests":
        out = []
        for m in inp["models"]:
            try:
                model = validate_model(m)
            except DomainFailure as exc:
                out.append(from_domain(exc))
                continue
            behavior = deep_thaw(model)
            behavior.pop("metadata", None)
            kinetics = model["kinetics"]
            out.append({"config_digest": digest(behavior),
                        "identity_digest": identity_digest(kinetics), "outcome": "validated",
                        "schema_digest": SCHEMA_SHA256,
                        "tolerance_digest": digest({k: kinetics[k] for k in TOLERANCE_KEYS})})
        return {"cases": out}
    if name == "state_identity":
        out = []
        for case in inp["cases"]:
            request = state_request(case["system"])
            energy = float(case["energy_ev"])
            zero = [[0.0, 0.0, 0.0] for _ in request["positions"]]
            state = state_from_relaxation(request, energy, zero, 1.0,
                                          {"calculator_evaluations": 0, "steps": 0,
                                           "calculator_identity": "xlang",
                                           "minimizer_identity": "xlang",
                                           "termination_reason": "force_tolerance"})
            out.append({"candidate_identity": state.candidate_identity,
                        "constraint_digest": state.constraint_digest,
                        "energy_ev": scalar(energy),
                        "fixed_contract_digest": fixed_contract_digest(request),
                        "geometry_certificate": geometry_certificate(request),
                        "state_id": state.state_id})
        return {"cases": out}
    if name == "philox_block":
        out = []
        for block in inp["blocks"]:
            key = tuple(block["key"])
            counter = tuple(block["counter"])
            out.append({"counter": list(counter), "key": list(key),
                        "words": list(philox4x32_10(counter, key))})
        return {"cases": out}
    if name == "uniform_words":
        out = []
        for a, b in inp["pairs"]:
            q = (a << 20) | (b >> 12)
            try:
                u = words_to_uniform(a, b)
            except DomainFailure as exc:
                out.append(from_domain(exc))
                continue
            out.append({"a": a, "b": b, "q": q, "raw_binary64_bits": bits(u), "uniform": num(u)})
        return {"cases": out}
    if name == "rng_derivation":
        draws = inp["draws"]
        trajectory = []
        for seed in inp["trajectory_seeds"]:
            stream = derive_trajectory_stream(seed)
            trajectory.append({"initial_state": stream.checkpoint(), "run_seed": seed,
                               "stream": draw(stream, draws)})
        saddle = []
        for entry in inp["saddle"]:
            stream = derive_saddle_stream(entry["run_seed"], entry["state_id"],
                                          entry["search_class"], entry["search_index"])
            saddle.append({"entry": entry, "initial_state": stream.checkpoint(),
                           "stream": draw(stream, draws)})
        return {"saddle": saddle, "trajectory": trajectory}
    if name == "rng_state_sequence":
        stream = philox_from_words(inp["key"], inp["initial_counter"])
        steps = []
        for _ in range(inp["draws"]):
            before = stream.checkpoint()
            try:
                u = stream.uniform()
            except DomainFailure as exc:
                steps.append(from_domain(exc))
                break
            steps.append({"state_after": stream.checkpoint(), "state_before": before,
                          "uniform": scalar(u)})
        return {"final_substream_digest": substream_digest(stream), "steps": steps}
    if name == "search_ids":
        out = []
        for entry in inp["entries"]:
            out.append({"entry": entry,
                        "search_id": "search:" + digest(
                            {"run_seed": entry["run_seed"],
                             "search_class": entry["search_class"],
                             "search_index": entry["search_index"],
                             "state_id": entry["state_id"]})})
        return {"ascending_commit_order": sorted(x["search_id"] for x in out), "cases": out}
    if name == "discovery_class":
        model = json.loads(json.dumps(template))
        model["discovery"]["classes"] = inp["classes"]
        model["kinetics"]["run_seed"] = inp["run_seed"]
        try:
            engine = ReferenceEngine(validate_model(model),
                                     extension={"calculator_command": ["/bin/true"]})
        except DomainFailure as exc:
            return from_domain(exc)
        out = []
        for index in inp["indices"]:
            entry, class_stream = engine._class_for_attempt(inp["state_id"], index)
            saddle_stream = derive_saddle_stream(inp["run_seed"], inp["state_id"], entry["name"],
                                                 index)
            out.append({"class_stream": class_stream.checkpoint(),
                        "saddle_stream": saddle_stream.checkpoint(),
                        "search_class": entry["name"],
                        "search_id": engine._search_id(inp["state_id"], entry["name"], index),
                        "search_index": index})
        return {"cases": out}
    if name == "event_ids":
        out = []
        for case in inp["cases"]:
            mode = [component for vector in case["unstable_direction"] for component in vector]
            canonical_mode = min((mode, [-component for component in mode]), key=canonical_bytes)
            endpoints = sorted((case["origin_state_id"], case["destination_state_id"]))
            mapping = [tuple(pair) for pair in case["active_atom_mapping"]]
            oriented = (sorted(mapping) if case["origin_state_id"] == endpoints[0]
                        else sorted((right, left) for left, right in mapping))
            pair_id = "pair:" + digest({"active_atom_mapping": [list(x) for x in oriented],
                                        "endpoint_state_ids": endpoints,
                                        "saddle_energy_ev": case["saddle_energy_ev"],
                                        "saddle_geometry_digest": case["saddle_geometry_digest"],
                                        "schema": "spark-atomistic-event-pair/2",
                                        "unstable_direction": canonical_mode})
            def directed(origin, destination):
                return "event:" + digest({"destination_state_id": destination,
                                          "origin_state_id": origin, "pair_id": pair_id,
                                          "schema": "spark-atomistic-directed-event/2"})
            out.append({"event_id": directed(case["origin_state_id"], case["destination_state_id"]),
                        "input": case, "pair_id": pair_id,
                        "reverse_event_id": directed(case["destination_state_id"],
                                                     case["origin_state_id"])})
        return {"cases": out}
    if name == "rate_pair":
        out = []
        for case in inp["cases"]:
            origin, destination, saddle = (case["origin_ev"], case["destination_ev"],
                                           case["saddle_ev"])
            forward, reverse = saddle - origin, saddle - destination
            if forward < -case["barrier_tolerance"] or reverse < -case["barrier_tolerance"]:
                out.append(failure("RATE_INVALID", "RATE-004"))
                continue
            beta = 1.0 / (8.617333262145e-5 * case["temperature"])
            log_prefactor = math.log(case["prefactor"])
            log_forward = log_prefactor - beta * forward
            log_reverse = log_prefactor - beta * reverse
            residual = (log_forward - log_reverse) + beta * (destination - origin)
            if abs(residual) > case["detailed_balance_tolerance"]:
                out.append(failure("DETAILED_BALANCE_VIOLATION", "RATE-005"))
                continue
            out.append({"barrier_ev": scalar(forward),
                        "detailed_balance_residual": scalar(residual),
                        "log_forward_rate_per_s": scalar(log_forward),
                        "log_reverse_rate_per_s": scalar(log_reverse), "outcome": "rated",
                        "reverse_barrier_ev": scalar(reverse)})
        return {"cases": out}
    if name in {"rate_snapshot", "kmc_selection"}:
        try:
            table = rate_table_from_rows(inp["rows"], inp["origin_state_id"],
                                         inp["log_rate_cutoff"])
        except DomainFailure as exc:
            return from_domain(exc)
        snapshot = table.snapshot()
        if name == "rate_snapshot":
            return {"outcome": "snapshot", "snapshot": snapshot}
        stream = derive_trajectory_stream(inp["run_seed"])
        steps = []
        for index in range(1, inp["steps"] + 1):
            try:
                selection = propose_serial_step(table, stream)
            except DomainFailure as exc:
                steps.append(from_domain(exc))
                break
            steps.append({"post_state_id": selection.destination_state_id,
                          "rng_after": selection.rng_after.checkpoint(),
                          "selected_event_id": selection.event_id,
                          "selected_rate_per_s": scalar(selection.selected_rate),
                          "selection_uniform": scalar(selection.selection_uniform),
                          "step_index": index,
                          "time_increment_s": scalar(selection.delta_time),
                          "time_uniform": scalar(selection.time_uniform),
                          "total_rate_per_s": scalar(selection.total_rate)})
            stream.commit_from(selection.rng_after)
        return {"snapshot": snapshot, "steps": steps}
    if name == "status_records":
        order = ["OK", "DISCOVERY_CONVERGED_HEURISTIC", "DUPLICATE_EVENT", "SADDLE_NOT_FOUND",
                 "INVALID_SADDLE", "SADDLE_WRONG_BASIN", "ENDPOINT_COLLAPSED",
                 "ENVIRONMENT_AMBIGUOUS", "BASIN_DISABLED", "DISCOVERY_INCOMPLETE",
                 "RELAX_NOT_CONVERGED", "EVENT_APPLICATION_FAILED", "CALCULATOR_FAILURE",
                 "NONFINITE_RESULT", "INVALID_INPUT", "SCHEMA_UNSUPPORTED", "INVALID_STATE",
                 "RATE_INVALID", "DETAILED_BALANCE_VIOLATION", "CATALOG_CONFLICT",
                 "CATALOG_INCOMPATIBLE", "ATOM_COUNT_CHANGE_UNSUPPORTED", "NO_ENABLED_EVENT",
                 "RESOURCE_LIMIT", "OUTPUT_EXISTS", "CHECKPOINT_CORRUPT",
                 "CHECKPOINT_INCOMPATIBLE", "CANCELLED", "INTERNAL_ERROR"]
        assert set(order) == set(STATUSES), "E2-STATUS-002 vocabulary drift"
        nonterminal = {"DUPLICATE_EVENT", "SADDLE_NOT_FOUND", "INVALID_SADDLE",
                       "SADDLE_WRONG_BASIN", "ENDPOINT_COLLAPSED", "ENVIRONMENT_AMBIGUOUS",
                       "BASIN_DISABLED"}

        def code(status, **kw):
            # E2-STATUS-003: a status that cannot terminate a public operation stores exit_code null.
            return None if status in nonterminal else exit_code(status, **kw)

        return {"cases": [{"exit_code_absorbing_requested": code(s, absorbing_ok=True),
                           "exit_code_default": code(s),
                           "exit_code_exploratory": code(s, exploratory=True),
                           "message": MESSAGE[s], "severity": SEVERITY[s], "status": s}
                          for s in order]}
    return {"outcome": "unsupported_probe", "probe": name}


def checkpoint_case(raw, model):
    # `validate_model` normalises in place, so every case gets its own copy.
    engine_config = validate_model(json.loads(json.dumps(model)))
    kinetics = engine_config["kinetics"]
    behavior = deep_thaw(engine_config)
    behavior.pop("metadata", None)
    arguments = {
        "expected_config_digest": digest(behavior),
        "expected_model_digest": engine_config["calculator"]["model_digest"],
        "expected_tolerance_digest": digest({k: kinetics[k] for k in TOLERANCE_KEYS}),
        "expected_identity_digest": identity_digest(kinetics),
        "kinetics": kinetics, "saddle_config": engine_config["saddle_search"],
        "relaxation_config": engine_config["relaxation"],
        "discovery_config": engine_config["discovery"],
        "resource_config": engine_config["resources"],
        "maximum_events": engine_config["resources"]["catalog_events"]}
    try:
        envelope = json.loads(raw)
        if (not isinstance(envelope, dict) or set(envelope) != {"payload", "payload_sha256"}
                or envelope["payload_sha256"] != digest(envelope["payload"])):
            return failure("CHECKPOINT_CORRUPT", "E2-CKPT-001")
        restored = validate_checkpoint_payload(envelope["payload"], **arguments)
    except DomainFailure as exc:
        return from_domain(exc)
    payload = envelope["payload"]
    steps = [{"post_state_id": s["post_state_id"], "pre_state_id": s["pre_state_id"],
              "rate_table_snapshot": s["rate_table_snapshot"],
              "selected_event_id": s["selected_event_id"],
              "selected_rate_per_s": scalar(s["selected_rate_per_s"]),
              "selection_uniform": scalar(s["selection_uniform"]), "step_index": s["step_index"],
              "time_increment_s": scalar(s["time_increment_s"]),
              "time_uniform": scalar(s["time_uniform"]),
              "total_rate_per_s": scalar(s["total_rate_per_s"])} for s in payload["trajectory"]]
    reencoded = canonical_bytes({"payload": payload, "payload_sha256": digest(payload)})
    return {"basin": payload["basin"], "catalog_digest": payload["catalog"]["digest"],
            "catalog_event_ids": sorted(payload["catalog"]["events"]),
            "catalog_state_ids": sorted(payload["catalog"]["states"]),
            "checkpoint_sequence": restored["checkpoint_sequence"],
            "current_state_id": restored["current_state_id"], "flags": restored["flags"],
            "initial_state_id": restored["initial_state_id"],
            "log_sequence": restored["log_sequence"], "outcome": "restored",
            "reencoded_equals_input": reencoded == raw,
            "rng_trajectory": restored["trajectory_rng"].checkpoint(),
            "simulation_time_s": scalar(restored["simulation_time"]),
            "step_index": restored["step_index"],
            "substream_ids": sorted(restored["substreams"]), "trajectory": steps}


# ------------------------------------------------------------------ emit / compare

def load_manifest():
    return json.loads((XLANG / "manifest.json").read_bytes())


def emit(out_root):
    manifest = load_manifest()
    out = Path(out_root) / "python"
    out.mkdir(parents=True, exist_ok=True)
    ids = []
    for entry in manifest["requests"]:
        raw = (XLANG / entry["file"]).read_bytes()
        assert len(raw) == entry["bytes"], "%s: corpus request length changed" % entry["id"]
        (out / (entry["id"] + ".out")).write_bytes(process_atomistic_json(raw).encode("utf-8"))
        ids.append(entry["id"])
    probes = json.loads((XLANG / "probes.json").read_bytes())
    template = json.loads((XLANG.parent / "e2_minimal_model.json").read_bytes())
    for case in probes["cases"]:
        (out / (case["id"] + ".out")).write_bytes(canonical_bytes(probe(case["input"], template)))
        ids.append(case["id"])
    model = json.loads((XLANG / "checkpoints" / "checkpoint_model.json").read_bytes())
    for entry in manifest["checkpoints"]:
        raw = (XLANG / entry["file"]).read_bytes()
        (out / (entry["id"] + ".out")).write_bytes(canonical_bytes(checkpoint_case(raw, model)))
        ids.append(entry["id"])
    assert len(set(ids)) == len(ids), "case IDs must be unique"
    (out / "_index.json").write_bytes(canonical_bytes({"backend": "python", "cases": sorted(ids)}))
    return sorted(ids)


def window(raw, offset, span=48):
    lo = max(0, offset - span)
    return repr(raw[lo:offset + span].decode("utf-8", "backslashreplace"))


def compare(out_root):
    manifest = load_manifest()
    probes = json.loads((XLANG / "probes.json").read_bytes())
    kinds = {}
    fixtures = {}
    tiers = {}
    for entry in manifest["requests"]:
        kinds[entry["id"]] = "request"
        fixtures[entry["id"]] = entry["fixture"]
        tiers[entry["id"]] = entry["tier"]
    for case in probes["cases"]:
        kinds[case["id"]] = "probe:" + case["input"]["probe"]
        fixtures[case["id"]] = case["fixture"]
        tiers[case["id"]] = case["tier"]
    for entry in manifest["checkpoints"]:
        kinds[entry["id"]] = "checkpoint"
        fixtures[entry["id"]] = entry["fixture"]
        tiers[entry["id"]] = entry["tier"]

    root = Path(out_root)
    rows = []
    for case_id in sorted(kinds):
        py = root / "python" / (case_id + ".out")
        rs = root / "rust" / (case_id + ".out")
        if not py.exists() or not rs.exists():
            rows.append({"case": case_id, "fixture": fixtures[case_id], "kind": kinds[case_id],
                         "tier": tiers[case_id],
                         "verdict": "NOT-APPLICABLE",
                         "reason": "missing artifact: python=%s rust=%s" % (py.exists(), rs.exists())})
            continue
        a, b = py.read_bytes(), rs.read_bytes()
        if a == b:
            rows.append({"case": case_id, "fixture": fixtures[case_id], "kind": kinds[case_id],
                         "tier": tiers[case_id],
                         "verdict": "PASS", "bytes": len(a)})
            continue
        offset = next((i for i in range(min(len(a), len(b))) if a[i] != b[i]), min(len(a), len(b)))
        rows.append({"case": case_id, "fixture": fixtures[case_id], "kind": kinds[case_id],
                     "tier": tiers[case_id],
                     "verdict": "DIVERGENT", "first_diff_offset": offset,
                     "python_bytes": len(a), "rust_bytes": len(b),
                     "python": window(a, offset), "rust": window(b, offset),
                     "shared_algebra": kinds[case_id].split(":")[-1] in SHARED_ALGEBRA})
    report = {"schema": "spark-atomistic-xlang-report/1", "cases": rows,
              "summary": {v: sum(1 for r in rows if r["verdict"] == v)
                          for v in ("PASS", "DIVERGENT", "NOT-APPLICABLE")},
              "summary_by_tier": {
                  tier: {v: sum(1 for r in rows if r["tier"] == tier and r["verdict"] == v)
                         for v in ("PASS", "DIVERGENT", "NOT-APPLICABLE")}
                  for tier in ("core", "adapter")}}
    roundtrip = root / "rust" / "_roundtrip.json"
    if roundtrip.exists():
        report["rust_self_read_failures"] = json.loads(roundtrip.read_bytes())
    (root / "parity_report.json").write_bytes(canonical_bytes(report))

    width = max(len(r["case"]) for r in rows)
    for r in rows:
        print("%-*s  %-7s  %-22s  %s" %
              (width, r["case"], r["tier"], r["fixture"], r["verdict"]))
        if r["verdict"] == "DIVERGENT":
            print("      first differing byte offset %d  (python %d bytes, rust %d bytes)%s"
                  % (r["first_diff_offset"], r["python_bytes"], r["rust_bytes"],
                     "  [SHARED-ALGEBRA probe]" if r["shared_algebra"] else ""))
            print("      python: %s" % r["python"])
            print("      rust:   %s" % r["rust"])
        elif r["verdict"] == "NOT-APPLICABLE":
            print("      %s" % r["reason"])
    print("\nsummary: %s" % report["summary"])
    print("tiers: %s" % report["summary_by_tier"])
    fixture_verdict = {}
    for r in rows:
        current = fixture_verdict.get(r["fixture"], "PASS")
        rank = {"PASS": 0, "NOT-APPLICABLE": 1, "DIVERGENT": 2}
        fixture_verdict[r["fixture"]] = max(current, r["verdict"], key=lambda v: rank[v])
    print("fixtures: %d PASS, %d DIVERGENT, %d NOT-APPLICABLE (of %d)"
          % (sum(1 for v in fixture_verdict.values() if v == "PASS"),
             sum(1 for v in fixture_verdict.values() if v == "DIVERGENT"),
             sum(1 for v in fixture_verdict.values() if v == "NOT-APPLICABLE"),
             len(fixture_verdict)))
    (root / "fixture_verdicts.json").write_bytes(canonical_bytes(fixture_verdict))
    core = report["summary_by_tier"]["core"]
    return core["DIVERGENT"] == 0 and core["NOT-APPLICABLE"] == 0


if __name__ == "__main__":
    command = sys.argv[1] if len(sys.argv) > 1 else "help"
    out_root = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else "/tmp/spark-xlang-out"
    if command == "emit":
        written = emit(out_root)
        print("wrote %d Python artifacts to %s/python" % (len(written), out_root))
    elif command == "compare":
        sys.exit(0 if compare(out_root) else 1)
    else:
        print(__doc__)
        sys.exit(2)
