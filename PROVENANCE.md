# SPARK source provenance

Status: package-scoped PolyForm Strict 1.0.0 adopted for the independent
atomistic Python/Rust candidates. Repository-wide release remains blocked by
rights review, contaminated Git history, and the missing root exact allowlist.

## Release boundary

The candidate core is the lattice KMC, dynamic-surface KMC, microkinetic,
reproducibility, statistics, and Python/Rust acceleration code developed in
this repository. It still requires file-level author and source review.
Published algorithms may guide behavior and mathematics; third-party source
code must not be copied or translated.

The former Python and Rust off-lattice modules were explicitly described as
ports of GPL-3.0 openFLY. They are excluded from the clean release boundary.
Renaming symbols, deleting comments, or translating a port does not make it an
independent implementation.

Independent Python and Rust replacements exist at
[`spark-atomistic`](spark-atomistic) and
[`spark-atomistic-rs`](spark-atomistic-rs). Their isolated implementers read
only the neutral specification with self-excluding SHA-256
`8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84`
and Errata 1 with SHA-256
`52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40`.
They do not read each other, the legacy off-lattice source, Git history,
quarantine, or competitor source. Errata 3 is adopted by D-127. Current gates:
Rust `cargo check --all-targets` has 0 errors/0 warnings and
`cargo test --all-targets` passes 69/69; Python unittest runs 100 tests with
99 passed and 1 cross-language comparison skipped when `SPARK_XLANG_OUT` is
absent; the separately executed cross-language gate passes 85/85 cases and
46/46 fixtures. No scientific validation or benchmark has run, so both remain
`implemented_unvalidated`, `validated=false`, `production=false`, and
`release=false`. They are not imported by the legacy SPARK root. Each package
contains the unmodified PolyForm Strict License 1.0.0 plus an explicit scope
notice; that grant does not cover repository history or excluded material.

The first independent Python candidate had one information-barrier exception:
a work-directory mistake exposed workspace file paths and separator-only
lines. It was retained as audit evidence under
`data/ai-kmc/quarantine/spark-atomistic-python-v1-cwd-incident-20260809` and is
permanently excluded from merge and release. The current Python candidate is a
fresh rewrite with zero reported barrier exceptions.

The first Rust candidate later reported a target-external `rg --files` path
exposure during static cleanup. No external file content was read or used, but
LIC-001 treats layout exposure as an exception. It is retained under
`data/ai-kmc/quarantine/spark-atomistic-rust-v1-layout-incident-20260809` and
permanently excluded. A fork-none Rust V2 implementer is rebuilding from the
same spec and erratum.

ReconKin and other repositories without an explicit license grant are used
only for public behavior and feature comparison. Their source is not a source
for SPARK implementation.

An independent graph-lattice candidate exists under
[`cleanroom/graph_kmc_v1`](cleanroom/graph_kmc_v1). Its implementer read only
the clean-room survey with full-file SHA-256
`37651ac4642df128e00cfde2e9a31ebe1a3bafd35c6b48933b578874dda19af1`
and the RNG erratum with full-file SHA-256
`fd6e5d5118b35ed9bd8702b19355ef573f7b569d18c0ae515a98a1b1bbc97002`.
It attests that no existing SPARK, competitor, quarantine, or Git-history
source was read. Its self-excluding file allowlist SHA-256 is
`b2cccf0b012e5a6b561626b55ae534cb4ae29b048b946a75caa02bded0291d78`.
The Python/Rust candidate remains isolated and `implemented_unvalidated`; no
import, Cargo build/check/test, golden run, or benchmark has executed.

The current multi-lattice implementation is not clean-room certified because
its design record says the implementer read GPL-3.0 kmcos source. It remains
outside a permissive release allowlist until rights or an independent
implementation record is established.

An isolated candidate now exists under
[`cleanroom/multilattice_v1`](cleanroom/multilattice_v1). Its implementer read
only the neutral paper specification with SHA-256
`db5fbe81317cbb25f5b8cf68bf6bd5e954a30d956980e30d15575dc293f61a23` and
attested that no existing SPARK or competitor implementation was read. It is
unintegrated and `implemented_unvalidated`. Isolated fix/review cycles reduced
static findings from `0 Critical / 5 High / 3 Medium`, then
`0 Critical / 1 High / 3 Medium`, to `0 Critical / 0 High / 0 Medium`.
Executed Python/Rust differential validation, contributor permission, and the
root license decision remain required.

The machine-readable quarantine is
[`provenance/QUARANTINE_V1.json`](provenance/QUARANTINE_V1.json). It also keeps
third-party benchmark outputs and PDFs as local evidence while excluding them
from a source release. Quarantine is conservative provenance control, not a
legal conclusion about any individual file.

## Contribution rule

Every new source file must be one of:

1. original work by a contributor who can grant the selected SPARK license;
2. a clean implementation from a written interface, paper, or standard, with
   the publication recorded in documentation; or
3. third-party code kept outside the SPARK release and governed by its own
   license.

Paper citation and acknowledgement remain required for scientific algorithms.
They do not replace source-license and copyright obligations.

Run `python3 tools/check_provenance.py` before building a release archive. It
checks the denylist, source-tree links/special files, candidate manifests,
hashes, and Rust clean-room assertions. It intentionally blocks public release
until an externally pinned, archive-enforced root exact allowlist exists. A
future passing result still cannot grant a license or establish copyright.
