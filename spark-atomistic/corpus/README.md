# Minimal corpus

`mock_calculator.py` is a deterministic stdlib-only calculator adapter used by the
conformance fixtures. Its analytic potential is `V = (x^2-1)^2 + y^2 + z^2` per atom in
eV with positions in Angstrom: minima at `x = +-1` (E = 0), a first-order saddle at
`x = 0` (E = 1.0 eV), and forward and reverse barriers both exactly 1.000 eV. It is a
conformance fixture, not a scientific benchmark.

Note that it raises `OverflowError` and exits nonzero once a coordinate grows past the
binary64 range. `ProcessCalculator` maps a nonzero exit to `CALCULATOR_FAILURE`
(`CALC-006`), which is a transaction-fail status rather than a candidate reject. See
`tests/TEST_SPEC.md` finding F3 for when that path is reachable.

The canonical model/request/response fixtures are NOT stored here. Errata 2
(`E2-PAR-001`) requires both backends to consume the same canonical corpus, so the
single copy lives at `../../spark-atomistic-rs/tests/corpus/` and is referenced from
`tests/fixtures.py` via the `CORPUS` constant and from `tests/xlang_harness.py` via
`XLANG`. Two shared-corpus problems measured on 2026-08-11 were reported in
`tests/TEST_SPEC.md` (findings F4 and F5) rather than edited there. The cross-language
sub-corpus `../../spark-atomistic-rs/tests/corpus/xlang/` was added on 2026-08-11 and is
described by its own README; it is the only shared-corpus material this work created,
and it exists because `E2-PAR-001` allows exactly one copy.

`minimal_model.json` was removed on 2026-08-11: it encoded the pre-Errata-2 layout
(`schema`/`system` envelope) and was rejected unconditionally by the current validator.
Because a rejected fixture still satisfies an `INVALID_INPUT` assertion, it silently
turned the retry/alpha rejection test into a false positive that passed without
exercising either rule. Fixtures asserted to be rejected MUST be paired with a baseline
that validates; every rejection test in `tests/` now carries that pairing.

The in-process `StubCalculator` in `tests/fixtures.py` reproduces this same potential
without the subprocess transport, so a solver experiment measures the solver.
