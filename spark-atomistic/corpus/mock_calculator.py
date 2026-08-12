# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
"""Deterministic analytic double-well JSON-lines calculator fixture."""

from __future__ import annotations

import hashlib
import json
import sys


MODEL_DIGEST = "7c799e3c0c25eb952d433430027d3d73de8d9f8f3d06064b6374f4b6eab4dd47"


def main() -> int:
    request = json.loads(sys.stdin.buffer.readline().decode("utf-8"))
    positions = request["state"]["positions"]
    energy = 0.0
    forces = []
    for x, y, z in positions:
        energy += (x * x - 1.0) ** 2 + y * y + z * z
        forces.append([-4.0 * x * (x * x - 1.0), -2.0 * y, -2.0 * z])
    response = {
        "status": "OK", "energy": energy, "forces": forces,
        "units": {"energy": "eV", "forces": "eV/angstrom"},
        "model_name": "analytic-double-well", "model_version": "1",
        "model_digest": MODEL_DIGEST,
        "evaluation_id": "eval-" + hashlib.sha256(
            request["request_digest"].encode("ascii")).hexdigest(),
        "deterministic": True, "request_digest": request["request_digest"],
    }
    sys.stdout.write(json.dumps(response, sort_keys=True, separators=(",", ":")) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

