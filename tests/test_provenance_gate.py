"""Minimal regression specification for the release provenance gate."""

from pathlib import Path

from tools.check_provenance import audit


def test_provenance_gate_accepts_original_core_and_rejects_port_markers(tmp_path):
    source = tmp_path / "spark"
    source.mkdir()
    clean = source / "core.py"
    clean.write_text("def rate():\n    return 1.0\n", encoding="utf-8")
    assert audit(Path(tmp_path)) == []

    clean.write_text("# Ported from openFLY\n", encoding="utf-8")
    assert audit(Path(tmp_path)) == [
        "forbidden provenance marker 'ported from openfly': spark/core.py"
    ]


def test_release_mode_retains_root_allowlist_blocker(tmp_path):
    assert audit(Path(tmp_path), release=True) == [
        "root release allowlist is not externally pinned or archive-enforced; public release blocked"
    ]
