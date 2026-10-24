"""Pure unit tests for multigene_identifier internals (blast stubbed)."""

from __future__ import annotations

import subprocess
import types
from pathlib import Path

import pytest

import gside.analysis.multigene_identifier as mg
from gside import config
from gside.analysis.multigene_identifier import (
    _blast_contigs,
    _db_version,
    _load_rules,
)


def _completed(stdout: str = "", returncode: int = 0, stderr: str = ""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


def _blastn_stub(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, script: str) -> Path:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    stub = bin_dir / "blastn"
    stub.write_text(script)
    stub.chmod(0o755)
    monkeypatch.setattr(config, "pixi_path", lambda: str(bin_dir))
    return stub


_HIT_LINE = "\t".join(
    [
        "q1",
        "mk~~~g1~~~x",
        "99.0",
        "400",
        "2",
        "0",
        "10",
        "409",
        "1",
        "400",
        "1e-50",
        "700",
        "500",
        "450",
    ]
)


class TestDbVersion:
    def test_missing_fasta_is_unknown(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mg, "_MARKERS_FASTA", tmp_path / "nope.fasta")
        assert _db_version() == "unknown"


class TestLoadRules:
    def test_missing_file_is_empty(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mg, "_MARKER_RULES", tmp_path / "nope.yaml")
        assert _load_rules() == []


class TestPriorityDeadCode:
    def test_priority_block_runs_without_effect(self, monkeypatch):
        monkeypatch.setattr(
            mg,
            "_load_rules",
            lambda: [
                {
                    "species": "X",
                    "genes": ["g1"],
                    "min_hits": 1,
                    "priority_over": ["Y"],
                }
            ],
        )
        monkeypatch.setattr(
            mg,
            "_blast_contigs",
            lambda contigs, threads=4: [{"gene": "g1", "identity": 99.0, "coverage": 100.0}],
        )
        from gside.analysis.multigene_identifier import identify_multigene

        assert identify_multigene("x.fna").to_dict()["species"] == "X"


class TestBlastParsing:
    def test_short_lines_skipped(self, monkeypatch, tmp_path):
        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\nexit 0\n")
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(f"short\n{_HIT_LINE}\n"))
        hits = _blast_contigs("x.fna")
        assert [h["gene"] for h in hits] == ["g1"]


class TestBlastFailureSurface:
    def test_nonzero_exit_raises_runtime_error(self, tmp_path, monkeypatch):
        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\necho 'BLAST engine error' >&2\nexit 3\n")
        with pytest.raises(RuntimeError, match=r"blastn failed \(exit 3\)"):
            _blast_contigs("x.fna")

    def test_failure_message_carries_stderr(self, tmp_path, monkeypatch):
        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\necho 'query file corrupt' >&2\nexit 2\n")
        with pytest.raises(RuntimeError, match="query file corrupt"):
            _blast_contigs("x.fna")

    def test_discovery_failure_raises_with_fix_hint(self, monkeypatch, tmp_path):
        monkeypatch.setattr(config, "pixi_path", lambda: str(tmp_path))
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("", 0))
        with pytest.raises(RuntimeError, match=r"blastn not found.*GSIDE_PIXI_BIN"):
            _blast_contigs("x.fna")

    def test_stub_hits_still_parsed(self, tmp_path, monkeypatch):
        _blastn_stub(
            tmp_path,
            monkeypatch,
            "#!/bin/sh\nprintf 'q1\\tmk~~~g1~~~x\\t99.0\\t400\\t2\\t0\\t10\\t409"
            "\\t1\\t400\\t1e-50\\t700\\t500\\t450\\n'\nexit 0\n",
        )
        hits = _blast_contigs("x.fna")
        assert [h["gene"] for h in hits] == ["g1"]


class TestBackendDelegation:
    def test_routes_through_blast_backend(self, tmp_path, monkeypatch):
        from gside.engine.backends import blast as blast_mod

        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\nexit 0\n")
        calls: dict[str, object] = {}

        class _B:
            def __init__(self, tool="blastn", threads=4, binary=None):
                calls["ctor"] = binary

            def find(self, query, db_path, **kw):
                calls["find"] = (query, db_path, kw)
                return []

        monkeypatch.setattr(blast_mod, "BlastBackend", _B)
        _blast_contigs("x.fna")
        assert calls["ctor"] == str(tmp_path / "bin" / "blastn")
        query, db_path, kw = calls["find"]
        assert query == Path("x.fna")
        assert str(db_path).endswith("markers_blastdb")
        assert kw["evalue"] == 1e-10
        assert kw["word_size"] == 11


class TestMainOutput:
    def test_json_and_text(self, monkeypatch, capsys, tmp_path):
        import gside.analysis.multigene_identifier as mg

        monkeypatch.setattr(
            mg,
            "identify_multigene",
            lambda c, threads=4: mg.MultiGeneResult(species="Salmonella", confidence="high"),
        )
        q = tmp_path / "q.fna"
        q.write_text(">q\nACGT\n")
        monkeypatch.setattr("sys.argv", ["prog", str(q), "--json"])
        mg.main()
        assert '"Salmonella"' in capsys.readouterr().out
        monkeypatch.setattr("sys.argv", ["prog", str(q)])
        mg.main()
        assert "Species: Salmonella (high)" in capsys.readouterr().out


class TestBestHitSelection:
    def test_short_fragment_does_not_displace_full_length_hit(self, monkeypatch):
        monkeypatch.setattr(
            mg,
            "_blast_contigs",
            lambda contigs, threads=4: [
                {
                    "gene": "ctxa",
                    "identity": 99.2,
                    "coverage": 12.5,
                    "ref_len": 2020,
                    "seqid": "markers~~~ctxa~~~x",
                },
                {
                    "gene": "ctxa",
                    "identity": 99.0,
                    "coverage": 100.0,
                    "ref_len": 2020,
                    "seqid": "markers~~~ctxa~~~x",
                },
            ],
        )
        result = mg.identify_multigene("q.fna")
        detected = {m["gene"]: m for m in result.detected_markers}
        assert "ctxa" in detected
        assert detected["ctxa"]["identity"] == 99.0


class TestRuleCoverageThresholds:
    def _identify(self, monkeypatch, rules, hits):
        monkeypatch.setattr(mg, "_load_rules", lambda: rules)
        monkeypatch.setattr(mg, "_blast_contigs", lambda contigs, threads=4: hits)
        from gside.analysis.multigene_identifier import identify_multigene

        return identify_multigene("x.fna").to_dict()

    def test_rule_min_coverage_below_global_is_honoured(self, monkeypatch):
        rules = [
            {
                "species": "X",
                "genes": ["g1"],
                "min_hits": 1,
                "min_identity": 90,
                "min_coverage": 50,
            }
        ]
        hits = [{"gene": "g1", "identity": 99.0, "coverage": 55.0}]
        assert self._identify(monkeypatch, rules, hits)["species"] == "X"

    def test_rule_without_min_coverage_uses_global_60(self, monkeypatch):
        rules = [{"species": "X", "genes": ["g1"], "min_hits": 1, "min_identity": 90}]
        hits = [{"gene": "g1", "identity": 99.0, "coverage": 55.0}]
        assert self._identify(monkeypatch, rules, hits)["species"] == "Unknown"

    def test_below_rule_min_coverage_is_rejected(self, monkeypatch):
        rules = [
            {
                "species": "X",
                "genes": ["g1"],
                "min_hits": 1,
                "min_identity": 90,
                "min_coverage": 50,
            }
        ]
        hits = [{"gene": "g1", "identity": 99.0, "coverage": 45.0}]
        assert self._identify(monkeypatch, rules, hits)["species"] == "Unknown"

    def test_lower_identity_hit_with_better_coverage_can_match(self, monkeypatch):
        rules = [
            {
                "species": "X",
                "genes": ["g1"],
                "min_hits": 1,
                "min_identity": 90,
                "min_coverage": 60,
            }
        ]
        hits = [
            {"gene": "g1", "identity": 99.0, "coverage": 40.0},
            {"gene": "g1", "identity": 95.0, "coverage": 90.0},
        ]
        assert self._identify(monkeypatch, rules, hits)["species"] == "X"


class TestGeneAlsoExcluded:
    def test_gene_in_both_lists_is_ignored(self, monkeypatch):
        monkeypatch.setattr(
            mg,
            "_load_rules",
            lambda: [
                {
                    "species": "X",
                    "genes": ["g1", "g2"],
                    "exclude_genes": ["g1"],
                    "min_hits": 1,
                }
            ],
        )
        monkeypatch.setattr(
            mg,
            "_blast_contigs",
            lambda contigs, threads=4: [{"gene": "g1", "identity": 99.0, "coverage": 100.0}],
        )
        from gside.analysis.multigene_identifier import identify_multigene

        assert identify_multigene("x.fna").to_dict()["species"] == "Unknown"


class TestFastqInput:
    def _rules(self, monkeypatch):
        monkeypatch.setattr(
            mg,
            "_load_rules",
            lambda: [{"species": "Salmonella", "genes": ["inva"], "min_hits": 1}],
        )

    def test_fastq_uses_minimap2_screening(self, monkeypatch, tmp_path):
        seen: dict[str, str] = {}

        def fake_screen(reads, markers, threads=4):
            seen["reads"] = str(reads)
            seen["markers"] = str(markers)
            return [{"gene": "inva", "identity": 99.0, "coverage": 100.0}]

        self._rules(monkeypatch)
        monkeypatch.setattr(mg, "screen_reads", fake_screen)
        fq = tmp_path / "reads.fq"
        fq.write_text("@r1\nACGT\n+\nIIII\n")
        from gside.analysis.multigene_identifier import identify_multigene

        assert identify_multigene(str(fq)).to_dict()["species"] == "Salmonella"
        assert seen["reads"] == str(fq)
        assert seen["markers"].endswith("markers.fasta")

    def test_fasta_still_uses_blast(self, monkeypatch, tmp_path):
        seen: dict[str, str] = {}

        def fake_blast(path, threads=4):
            seen["path"] = str(path)
            return [{"gene": "inva", "identity": 99.0, "coverage": 100.0}]

        self._rules(monkeypatch)
        monkeypatch.setattr(mg, "_blast_contigs", fake_blast)
        genome = tmp_path / "g.fna"
        genome.write_text(">c1\nACGT\n")
        from gside.analysis.multigene_identifier import identify_multigene

        assert identify_multigene(str(genome)).to_dict()["species"] == "Salmonella"
        assert seen["path"] == str(genome)


def _pick(*genes: str) -> str:
    """Species chosen by the real rules for a synthetic set of detected genes."""
    rules = _load_rules()
    by_gene = {gene: [{"gene": gene, "identity": 99.0, "coverage": 100.0}] for gene in genes}
    significant = {gene: by_gene[gene][0] for gene in genes}
    species, _confidence, _rule = mg._pick_rule(
        rules, by_gene, significant, mg._priority_species(rules)
    )
    return species


class TestRuleDesign:
    def test_shigella_genus_is_uida_ipah_lacy_negative(self):
        # No species marker -> genus-level Shigella (not forced to boydii).
        assert _pick("uida", "ipah") == "Shigella"
        assert _pick("uida", "ipah", "lacy") == "Escherichia_coli"

    def test_shigella_flexneri_uses_rfc_or_rfb_marker(self):
        assert _pick("uida", "ipah", "rfc") == "Shigella_flexneri"
        assert _pick("uida", "ipah", "rfbe") == "Shigella_flexneri"

    def test_shigella_sonnei_uses_wbgz_or_chromosomal_marker(self):
        assert _pick("uida", "ipah", "wbgz") == "Shigella_sonnei"
        assert _pick("uida", "ipah", "sonnei_psuk") == "Shigella_sonnei"

    def test_shigella_dysenteriae_uses_rfpb_or_sd_panel(self):
        assert _pick("uida", "ipah", "rfpb") == "Shigella_dysenteriae"
        assert _pick("uida", "ipah", "sd1_wzx") == "Shigella_dysenteriae"

    def test_shigella_boydii_uses_sb_serotype_panel(self):
        assert _pick("uida", "ipah", "sb2_wzx") == "Shigella_boydii"
        assert _pick("uida", "ipah", "rfc") == "Shigella_flexneri"

    def test_eiec_is_escherichia_coli_not_shigella(self):
        assert _pick("uida", "ipah", "lacy") == "Escherichia_coli"

    def test_bordetella_pertussis_needs_both_targets(self):
        assert _pick("ptxs1", "tcfa") == "Bordetella_pertussis"
        assert _pick("ptxs1") == "Unknown"

    def test_renamed_labels(self):
        assert _pick("tlh") == "Vibrio_parahaemolyticus"
        assert _pick("inva") == "Salmonella_enterica"


def test_best_per_gene_prefers_coverage_on_identity_tie():
    hits = [
        {"gene": "ipah", "identity": 100.0, "coverage": 74.1, "ref_len": 100},
        {"gene": "ipah", "identity": 100.0, "coverage": 100.0, "ref_len": 100},
    ]
    assert mg._best_per_gene(hits)["ipah"]["coverage"] == 100.0


def test_detected_marker_ties_are_gene_sorted(monkeypatch, tmp_path):
    monkeypatch.setattr(mg, "_load_rules", lambda: [])
    monkeypatch.setattr(mg, "_db_version", lambda: "test")
    monkeypatch.setattr(
        mg, "_blast_contigs", lambda path, threads: [
            {"gene": "b", "identity": 100.0, "coverage": 100.0, "ref_len": 100},
            {"gene": "a", "identity": 100.0, "coverage": 100.0, "ref_len": 100},
        ]
    )
    # No rule is required to exercise detected_markers ordering.
    result = mg.identify_multigene(str(tmp_path / "g.fna"))
    assert [m["gene"] for m in result.detected_markers] == ["a", "b"]
