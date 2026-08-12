#!/usr/bin/env python3
"""Audit SPARK clean-room integrity; optionally enforce release-only gates."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import unicodedata


SOURCE_SUFFIXES = frozenset({".py", ".rs"})
SOURCE_ROOTS = (
    "spark",
    "spark-rs/src",
    "examples",
    "cleanroom",
    "spark-atomistic",
    "spark-atomistic-rs",
)
FORBIDDEN_PATHS = (
    "spark/offlattice",
    "spark-rs/src/offlattice",
    "spark-rs/src/python_bindings.rs",
    "examples/offlattice_fe_vacancy.py",
    "docs/reference-software-summary/SOURCE_CODE_ANALYSIS.md",
)
FORBIDDEN_MARKERS = (
    "ported from openfly",
    "port of openfly",
    "matches openfly exactly",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _valid_digest(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value)
    )


def _valid_manifest_path(value: object) -> bool:
    if (
        not isinstance(value, str)
        or not value
        or "\\" in value
        or unicodedata.normalize("NFC", value) != value
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        return False
    path = PurePosixPath(value)
    return (
        not path.is_absolute()
        and bool(path.parts)
        and value == path.as_posix()
        and all(part not in ("", ".", "..") for part in path.parts)
    )


def _portable_collisions(paths: set[str], label: str) -> list[str]:
    seen: dict[str, str] = {}
    findings: list[str] = []
    for path in sorted(paths):
        portable = unicodedata.normalize("NFC", path).casefold()
        previous = seen.get(portable)
        if previous is not None and previous != path:
            findings.append(f"portable-name collision in {label}: {previous!r} / {path!r}")
        seen[portable] = path
    return findings


def _strict_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _symlink_findings(package: Path, label: str) -> list[str]:
    findings: list[str] = []
    for path in sorted(package.rglob("*")):
        relative = path.relative_to(package).as_posix()
        if path.is_symlink():
            findings.append(f"symlink forbidden in {label}: {relative}")
        elif not path.is_dir() and not path.is_file():
            findings.append(f"special file forbidden in {label}: {relative}")
    return findings


def _audit_python_manifest(root: Path) -> list[str]:
    package = root / "spark-atomistic"
    manifest = package / "FILE_HASHES.sha256"
    if not package.is_dir() or package.is_symlink():
        return ["required atomistic Python package missing or unsafe: spark-atomistic"]
    if not manifest.is_file() or manifest.is_symlink():
        return ["missing atomistic Python hash manifest: spark-atomistic/FILE_HASHES.sha256"]

    expected: dict[str, str] = {}
    findings = _symlink_findings(package, "atomistic Python package")
    for line_number, raw in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        parts = raw.split(maxsplit=1)
        if len(parts) != 2 or not _valid_digest(parts[0]):
            findings.append(f"invalid Python hash manifest line {line_number}")
            continue
        digest, relative = parts[0].lower(), parts[1].strip()
        if not _valid_manifest_path(relative):
            findings.append(f"invalid Python hash manifest path: {relative!r}")
            continue
        if relative in expected:
            findings.append(f"duplicate Python hash manifest path: {relative}")
        expected[relative] = digest

    actual = {
        path.relative_to(package).as_posix()
        for path in package.rglob("*")
        if (
            path.is_file()
            and not path.is_symlink()
            and path != manifest
            and "__pycache__" not in path.parts
            and path.suffix != ".pyc"
        )
    }
    findings.extend(_portable_collisions(set(expected), "Python manifest"))
    findings.extend(_portable_collisions(actual, "atomistic Python package"))
    for relative in sorted(actual - expected.keys()):
        findings.append(f"unlisted atomistic Python file: spark-atomistic/{relative}")
    for relative in sorted(expected.keys() - actual):
        findings.append(f"missing atomistic Python file: spark-atomistic/{relative}")
    for relative in sorted(actual & expected.keys()):
        if _sha256(package / relative) != expected[relative]:
            findings.append(f"atomistic Python hash mismatch: spark-atomistic/{relative}")
    return findings


def _audit_rust_manifest(root: Path) -> list[str]:
    package = root / "spark-atomistic-rs"
    manifest = package / "PROVENANCE.json"
    if not package.is_dir() or package.is_symlink():
        return ["required atomistic Rust package missing or unsafe: spark-atomistic-rs"]
    if not manifest.is_file() or manifest.is_symlink():
        return ["missing atomistic Rust provenance manifest: spark-atomistic-rs/PROVENANCE.json"]

    try:
        data = json.loads(
            manifest.read_text(encoding="utf-8"), object_pairs_hook=_strict_object
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        return [f"invalid atomistic Rust provenance manifest: {error}"]
    if not isinstance(data, dict):
        return ["invalid atomistic Rust provenance manifest: root must be object"]

    expected: dict[str, str] = {}
    findings = _symlink_findings(package, "atomistic Rust package")
    required_top_level = {
        "schema",
        "package",
        "backend",
        "ir",
        "claims",
        "allowed_sources",
        "attestation",
        "dynamic_actions",
        "dependencies",
        "project_license",
        "project_license_assigned",
        "review_repairs",
        "sha256_allowlist",
        "self_hash_rule",
        "cross_language_byte_identity",
        "errata_3_candidate",
    }
    if set(data) != required_top_level:
        findings.append("invalid Rust provenance top-level key set")
    if data.get("schema") != "spark-clean-room-provenance/1":
        findings.append("invalid Rust provenance schema")
    if (
        data.get("package") != "spark-atomistic-rs"
        or data.get("backend") != "rust"
        or data.get("ir") != "spark-atomistic-model/1"
    ):
        findings.append("invalid Rust package/backend/IR provenance")
    required_claims = {
        "validated": False,
        "production": False,
        "release": False,
    }
    if data.get("claims") != required_claims:
        findings.append("invalid or overclaimed Rust implementation provenance")
    expected_sources = [
        {
            "path": "/home/shidi/ai-chemist/catgo-projects/data/ai-kmc/specs/OFFLATTICE_OTF_KMC_SPEC_V1.md",
            "sha256": "8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84",
            "hash_scope": "normative self-excluding digest stated by source",
        },
        {
            "path": "/home/shidi/ai-chemist/catgo-projects/data/ai-kmc/specs/OFFLATTICE_OTF_KMC_SPEC_V1_ERRATA_1.md",
            "sha256": "52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40",
            "hash_scope": "normative self-excluding digest stated by source",
        },
        {
            "path": "/home/shidi/ai-chemist/catgo-projects/data/ai-kmc/specs/OFFLATTICE_OTF_KMC_SPEC_V1_ERRATA_2_PARITY.md",
            "sha256": "eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995",
            "hash_scope": "normative self-excluding digest stated by source",
        },
        {
            "path": "/home/shidi/ai-chemist/catgo-projects/data/ai-kmc/specs/OFFLATTICE_OTF_KMC_SPEC_V1_ERRATA_3.md",
            "sha256": "eba384af3694c5f3997caf28829e56d188ef9929f29d99a2116520f0067d8a96",
            "hash_scope": "normative self-excluding digest stated by source",
        },
    ]
    if data.get("allowed_sources") != expected_sources:
        findings.append("unexpected Rust implementer source provenance")
    expected_attestation = {
        "independently_authored": True,
        "prohibited_sources_read": 0,
        "source_code_sources_consulted": 0,
        "tests_copied": 0,
        "symbols_copied": 0,
        "file_layout_copied": 0,
        "target_external_metadata_exposed": False,
    }
    if data.get("attestation") != expected_attestation:
        findings.append("invalid Rust clean-room attestation")
    dynamic_actions = data.get("dynamic_actions")
    if not isinstance(dynamic_actions, dict) or set(dynamic_actions) != {
        "cargo_build",
        "cargo_check",
        "cargo_test",
        "rustfmt",
        "implementation_executed",
        "benchmark_run",
    } or any(not isinstance(value, bool) for value in dynamic_actions.values()) or dynamic_actions.get("benchmark_run") is not False:
        findings.append("invalid Rust dynamic-action provenance")
    expected_license = {
        "decision": "D-130",
        "license_file": "LICENSE",
        "license_sha256": "e2361f52ad5be22b937a6e983c824a534c5cffa454b6c34af2f8ce0c2cdf7c1a",
        "name": "PolyForm Strict License 1.0.0",
        "scope_file": "LICENSE_SCOPE.md",
        "source": "https://polyformproject.org/licenses/strict/1.0.0.txt",
    }
    if data.get("project_license_assigned") is not True or data.get("project_license") != expected_license:
        findings.append("invalid Rust project-license provenance")
    expected_dependencies = [
        {
            "name": "serde",
            "version": "1.0.219",
            "declared_license": "MIT OR Apache-2.0",
            "linking": "Rust static dependency",
            "audited": False,
        },
        {
            "name": "serde_json",
            "version": "1.0.140",
            "declared_license": "MIT OR Apache-2.0",
            "linking": "Rust static dependency",
            "audited": False,
        },
        {
            "name": "sha2",
            "version": "0.10.8",
            "declared_license": "MIT OR Apache-2.0",
            "linking": "Rust static dependency",
            "audited": False,
        },
    ]
    if data.get("dependencies") != expected_dependencies:
        findings.append("invalid Rust dependency provenance")
    expected_repairs = [
        *(f"H{index}" for index in range(1, 13)),
        *(f"M{index}" for index in range(1, 8)),
        *(f"R2-H{index}" for index in range(1, 6)),
        *(f"R2-M{index}" for index in range(1, 8)),
        "R3-H1",
        *(f"R3-M{index}" for index in range(1, 6)),
        "R4-M1",
        "R4-M2",
    ]
    repairs = data.get("review_repairs")
    if (
        not isinstance(repairs, list)
        or len(repairs) != len(set(repairs))
        or any(not isinstance(item, str) or not item for item in repairs)
        or not set(expected_repairs).issubset(repairs)
    ):
        findings.append("invalid Rust review-repair provenance")
    if data.get("self_hash_rule") != (
        "PROVENANCE.json is excluded from its own allowlist; every other "
        "package file is mandatory and no unlisted file is allowed."
    ):
        findings.append("invalid Rust self-hash rule")
    values = data.get("sha256_allowlist")
    if not isinstance(values, dict):
        findings.append("missing Rust provenance allowlist: sha256_allowlist")
        values = {}
    for relative, digest in values.items():
        if not _valid_manifest_path(relative):
            findings.append(f"invalid Rust provenance path: {relative!r}")
            continue
        if relative in expected:
            findings.append(f"duplicate Rust provenance path: {relative}")
        if not _valid_digest(digest):
            findings.append(f"invalid Rust provenance digest: {relative}")
            continue
        expected[relative] = digest.lower()

    actual = {
        path.relative_to(package).as_posix()
        for path in package.rglob("*")
        if path.is_file() and not path.is_symlink() and path != manifest
    }
    findings.extend(_portable_collisions(set(expected), "Rust provenance manifest"))
    findings.extend(_portable_collisions(actual, "atomistic Rust package"))
    for relative in sorted(actual - expected.keys()):
        findings.append(f"unlisted atomistic Rust file: spark-atomistic-rs/{relative}")
    for relative in sorted(expected.keys() - actual):
        findings.append(f"missing atomistic Rust file: spark-atomistic-rs/{relative}")
    for relative in sorted(actual & expected.keys()):
        if _sha256(package / relative) != expected[relative]:
            findings.append(f"atomistic Rust hash mismatch: spark-atomistic-rs/{relative}")
    return findings


def audit(root: Path, *, release: bool = False) -> list[str]:
    root = root.resolve()
    findings: list[str] = []

    for relative in FORBIDDEN_PATHS:
        candidate = root / relative
        if candidate.exists() or candidate.is_symlink():
            findings.append(f"forbidden release path exists: {relative}")

    if release:
        findings.append(
            "root release allowlist is not externally pinned or archive-enforced; "
            "public release blocked"
        )

    for source_root in SOURCE_ROOTS:
        base = root / source_root
        if not base.exists() and not base.is_symlink():
            continue
        if base.is_symlink() or not base.is_dir():
            findings.append(f"unsafe source root: {source_root}")
            continue
        for path in sorted(base.rglob("*")):
            relative = path.relative_to(root).as_posix()
            if path.is_symlink():
                findings.append(f"symlink forbidden in source tree: {relative}")
                continue
            if not path.is_dir() and not path.is_file():
                findings.append(f"special file forbidden in source tree: {relative}")
                continue
            if not path.is_file() or path.suffix not in SOURCE_SUFFIXES:
                continue
            text = path.read_text(encoding="utf-8", errors="replace").lower()
            for marker in FORBIDDEN_MARKERS:
                if marker in text:
                    findings.append(
                        f"forbidden provenance marker {marker!r}: {relative}"
                    )
    if (root / "spark-atomistic").exists():
        findings.extend(_audit_python_manifest(root))
    if (root / "spark-atomistic-rs").exists():
        findings.extend(_audit_rust_manifest(root))
    return findings


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    release = sys.argv[1:] == ["--release"]
    if sys.argv[1:] not in ([], ["--release"]):
        print("usage: check_provenance.py [--release]", file=sys.stderr)
        return 2
    findings = audit(root, release=release)
    if findings:
        for finding in findings:
            print(f"PROVENANCE_FAIL: {finding}", file=sys.stderr)
        return 1
    print("PROVENANCE_CLEANROOM_MANIFESTS_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
