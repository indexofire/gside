"""Pure unit tests for engine registry, backend registry, read sniffing."""

from __future__ import annotations

import gzip

import pytest

from gside.engine.backends import available, get_backend, register
from gside.engine.read_mapper import ReadMapper, _sniff_read_type
from gside.engine.registry import Registry


class TestRegistry:
    def test_case_insensitive(self):
        reg = Registry()
        reg.register("BlastN", object)
        assert reg.get("blastn") is object
        assert reg.has(" BLASTN ")

    def test_empty_name(self):
        with pytest.raises(ValueError):
            Registry().register("  ", object)

    def test_unknown(self):
        with pytest.raises(KeyError):
            Registry().get("nope")

    def test_available_snapshot(self):
        reg = Registry()
        reg.register("a", object)
        assert reg.available() == {"a": object}


class TestBackendRegistry:
    def test_available_lists_all_builtins(self):
        names = available()
        assert {"blastn", "skani", "mash", "sourmash", "minimap2", "mmseqs2"} <= set(names)

    def test_unknown_backend(self):
        with pytest.raises(KeyError):
            get_backend("nope-not-real")

    def test_kma_backend_resolves(self, monkeypatch):
        import gside.engine.backends.kma as kma_mod
        from gside.engine.backends.kma import KmaBackend

        monkeypatch.setattr(kma_mod, "which", lambda tool: "/bin/kma")
        assert isinstance(get_backend("kma"), KmaBackend)

    def test_custom_register(self):
        class _Dummy:
            pass

        register("unit-dummy-xyz", _Dummy)
        assert isinstance(get_backend("unit-dummy-xyz"), _Dummy)


def _fastq(path, seq_len: int, gz: bool = False) -> str:
    seq = "A" * seq_len
    text = f"@r1\n{seq}\n+\n{'I' * seq_len}\n"
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return path


class TestSniff:
    def test_empty_list(self):
        assert _sniff_read_type([]) == "short"

    def test_missing_file(self, tmp_path):
        assert _sniff_read_type([str(tmp_path / "nope.fq")]) == "short"

    def test_short_reads(self, tmp_path):
        assert _sniff_read_type([_fastq(str(tmp_path / "a.fq"), 150)]) == "short"

    def test_long_reads(self, tmp_path):
        assert _sniff_read_type([_fastq(str(tmp_path / "b.fq"), 2000)]) == "long"

    def test_gzipped(self, tmp_path):
        assert _sniff_read_type([_fastq(str(tmp_path / "c.fq.gz"), 2000, gz=True)]) == "long"


class TestSelect:
    def test_fasta_goes_minimap2(self):
        assert ReadMapper._select(["x.fna"]) == "minimap2"

    def test_explicit_type(self):
        assert ReadMapper._select(["x.fq"], read_type="long") == "minimap2"
        assert ReadMapper._select(["x.fq"], read_type="short") == "bwa"

    def test_sniffed(self, tmp_path):
        assert ReadMapper._select([_fastq(str(tmp_path / "d.fq"), 150)]) == "bwa"


class TestBlastpToolKwarg:
    def test_tool_injected(self):
        import gside.engine.backends as backends

        seen = {}

        class _Dummy:
            def __init__(self, **kw):
                seen.update(kw)

        backends.register("blastp", _Dummy)
        try:
            backends.get_backend("blastp")
        finally:
            del backends._REG._store["blastp"]
        assert seen.get("tool") == "blastp"
