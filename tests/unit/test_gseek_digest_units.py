"""Real _reference_digest branches (file vs directory reference).

Kept in its own module because test_gseek_bridge_units autouse-patches
``_reference_digest`` for every test there; these tests exercise the real
implementation against tmp_path files and directories.
"""

from __future__ import annotations

import pytest

from gside.engine import gseek_bridge as bridge


@pytest.fixture(autouse=True)
def _clear_digest_cache():
    bridge._REFERENCE_DIGESTS.clear()
    yield
    bridge._REFERENCE_DIGESTS.clear()


def test_reference_digest_of_a_file_uses_size_and_sha256(tmp_path):
    ref = tmp_path / "ref.msh"
    ref.write_bytes(b"abc")
    size, digest = bridge._reference_digest(ref)
    assert size == 3
    assert len(digest) == 64


def test_reference_digest_of_a_directory_uses_tree_size_and_tree_hash(tmp_path):
    ref_dir = tmp_path / "panel.sketch"
    ref_dir.mkdir()
    (ref_dir / "a.sketch").write_bytes(b"xyz")
    (ref_dir / "b.sketch").write_bytes(b"wxyz")
    size, digest = bridge._reference_digest(ref_dir)
    assert size == 7
    assert len(digest) == 64


def test_reference_digest_is_cached_per_resolved_path(tmp_path):
    ref = tmp_path / "ref.msh"
    ref.write_bytes(b"abc")
    first = bridge._reference_digest(ref)
    ref.write_bytes(b"different content")
    second = bridge._reference_digest(ref)
    assert first == second
