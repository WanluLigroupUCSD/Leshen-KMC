# Off-lattice/on-the-fly KMC specification v1 — Errata 1

Date: 2026-08-09  
Status: normative  
Applies to: `OFFLATTICE_OTF_KMC_SPEC_V1.md`  
Base self-excluding SHA-256: `8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84`

## Correction

The original `DET-002` mapping can round to binary64 `1.0` when `m=2^53-1`. It is replaced by the requirements below.

`E1-DET-002-A` Each binary64 uniform MUST consume one ordered pair of unsigned 32-bit output words from `Philox4x32-10`. For each output block, words are ordered by algorithm lane as `w0, w1, w2, w3`; the first uniform uses `(w0,w1)` and the second uses `(w2,w3)`. Blocks remain in increasing counter order. Host byte order MUST NOT alter lane or bit order.

`E1-DET-002-B` For an ordered pair `(a,b)`, construct the unsigned 52-bit integer

\[
q=(a\ll20)\;|\;(b\gg12),\qquad 0\le q<2^{52}.
\]

Thus `a[31:0]` becomes `q[51:20]`, `b[31:12]` becomes `q[19:0]`, and `b[11:0]` is discarded.

`E1-DET-002-C` The returned value MUST be the exactly representable binary64 midpoint

\[
u=(2q+1)2^{-53}.
\]

An equivalent exact bit construction MAY be used only when it produces the identical IEEE-754 binary64 bit pattern. Consequently `0<u<1` for every output.

`E1-DET-002-D` The boundary golden values are normative:

| `q` | Exact value | Hexadecimal binary64 | Raw binary64 bits |
|---:|---|---|---|
| `0` | `2^-53` | `0x1.0000000000000p-53` | `0x3ca0000000000000` |
| `2^52-1` | `1-2^-53` | `0x1.fffffffffffffp-1` | `0x3fefffffffffffff` |

`E1-DET-002-E` Python and Rust golden fixtures MUST use this word pairing, bit selection, mapping, and boundary table. The generated uniform bit patterns MUST be byte-identical. All other `PAR-*` tolerances and requirements remain unchanged.

## Effect

This erratum normatively overrides `DET-002` and every parity/golden-fixture expectation derived from it. The base specification file and its self-excluding SHA-256 remain unchanged. No other requirement, interface, status, priority, or scope changes.

## Erratum digest

Hash rule: SHA-256 of the exact UTF-8 bytes before the line beginning `## Erratum digest`; the separator blank line is included and this digest section is excluded.

Erratum SHA-256: `52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40`
