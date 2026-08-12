# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""One canonical request on stdin, one canonical response on stdout."""

from __future__ import annotations

import json
import sys

from .api import process_atomistic_json


def main() -> int:
    response = process_atomistic_json(sys.stdin.buffer.read())
    sys.stdout.write(response)
    try:
        decoded = json.loads(response)
        code = decoded.get("exit_code", 70)
        return code if type(code) is int else 70
    except (json.JSONDecodeError, AttributeError):
        return 70


if __name__ == "__main__":
    raise SystemExit(main())
