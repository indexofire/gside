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

    def test_custom_marker_db_builds_its_own_blast_db(self, tmp_path, monkeypatch):
        from gside.engine.backends import blast as blast_mod
        from gside.marker_db import MarkerDatabase

        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\nexit 0\n")
        root = tmp_path / "mydb"
        root.mkdir()
        (root / "markers.fasta").write_text(">markers~~~g1~~~x marker\nACGTACGT\n")
        (root / "marker_rules.yaml").write_text(
            "rules:\n- species: Test_species\n  genes: [g1]\n  min_hits: 1\n"
        )
        calls: dict[str, object] = {}

        class _B:
            def __init__(self, tool="blastn", threads=4, binary=None):
                calls["ctor"] = binary

            def make_db(self, fasta, prefix):
                calls["make_db"] = (fasta, prefix)
                Path(f"{prefix}.nhr").write_bytes(b"")
                return prefix

            def find(self, query, db_path, **kw):
                calls["find_db"] = str(db_path)
                return []

        monkeypatch.setattr(blast_mod, "BlastBackend", _B)
        db = MarkerDatabase.from_path(root)
        _blast_contigs("x.fna", 4, db)
        fasta, prefix = calls["make_db"]
        assert fasta == root / "markers.fasta"
        assert str(prefix).startswith(str(root))
        assert calls["find_db"].startswith(str(root))

    def test_shadow_engine_runs_both_and_compares(self, tmp_path, monkeypatch):
        import types as _types

        from gside.engine import gseek_bridge

        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\nexit 0\n")
        monkeypatch.setenv("GSIDE_ENGINE", "shadow")

        def hit(subject):
            return _types.SimpleNamespace(
                subject_id=subject, identity=99.0, subject_coverage=100.0, subject_length=400
            )

        monkeypatch.setattr(mg, "_legacy_blast_hits", lambda c, t, db=None: [hit("mk~~~g1~~~x")])
        monkeypatch.setattr(gseek_bridge, "blast_hits", lambda q, m, **kw: [hit("mk~~~g1~~~x")])
        seen: dict[str, object] = {}
        monkeypatch.setattr(
            gseek_bridge, "log_blast_shadow", lambda old, new: seen.update(n=len(new))
        )
        hits = _blast_contigs("x.fna")
        assert [h["gene"] for h in hits] == ["g1"]
        assert seen == {"n": 1}


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

        def fake_screen(reads, markers, threads=4, min_mapq=20):
            seen["reads"] = str(reads)
            seen["markers"] = str(markers)
            seen["min_mapq"] = min_mapq
            return [{"gene": "inva", "identity": 99.0, "coverage": 100.0}]

        self._rules(monkeypatch)
        monkeypatch.setattr(mg, "screen_reads", fake_screen)
        fq = tmp_path / "reads.fq"
        fq.write_text("@r1\nACGT\n+\nIIII\n")
        from gside.analysis.multigene_identifier import identify_multigene

        assert identify_multigene(str(fq)).to_dict()["species"] == "Salmonella"
        assert seen["reads"] == str(fq)
        assert seen["markers"].endswith("markers.fasta")
        assert seen["min_mapq"] == 20

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

    def test_reads_mode_reports_control_and_gates_low_depth(self, monkeypatch, tmp_path):
        self._rules(monkeypatch)
        monkeypatch.setattr(mg, "_control_config", lambda db=None: ({"Salmonella": ["rpod"]}, {}))
        monkeypatch.setattr(
            mg,
            "screen_reads",
            lambda reads, markers, threads=4, min_mapq=20: [
                {"gene": "inva", "identity": 99.0, "coverage": 100.0, "mean_depth": 30.0},
                {"gene": "hilA", "identity": 99.0, "coverage": 100.0, "mean_depth": 2.0},
                {"gene": "rpod", "identity": 99.0, "coverage": 100.0, "mean_depth": 40.0},
            ],
        )
        fq = tmp_path / "reads.fq"
        fq.write_text("@r1\nACGT\n+\nIIII\n")
        from gside.analysis.multigene_identifier import identify_multigene

        payload = identify_multigene(str(fq)).to_dict()
        assert payload["species"] == "Salmonella"
        assert payload["control_gene"] == "rpod"
        assert payload["control_depth"] == 40.0
        assert payload["control_source"] == "dedicated"
        assert payload["depth_gate"] == (
            "min_depth=5x; min_relative_depth=10%; min_mapq=20; control=rpod"
        )
        markers = {m["gene"]: m for m in payload["detected_markers"]}
        assert set(markers) == {"inva"}  # hilA gated out, rpod is a hidden control
        assert markers["inva"]["mean_depth"] == 30.0
        assert markers["inva"]["relative_depth_pct"] == 75.0
        assert any("hilA" in note for note in payload["notes"])

    def test_reads_mode_without_control_uses_sample_median(self, monkeypatch, tmp_path):
        self._rules(monkeypatch)
        monkeypatch.setattr(mg, "_control_config", lambda db=None: ({}, {}))
        monkeypatch.setattr(
            mg,
            "screen_reads",
            lambda reads, markers, threads=4, min_mapq=20: [
                {"gene": "inva", "identity": 99.0, "coverage": 100.0, "mean_depth": 10.0},
                {"gene": "hilA", "identity": 99.0, "coverage": 100.0, "mean_depth": 30.0},
            ],
        )
        fq = tmp_path / "reads.fq"
        fq.write_text("@r1\nACGT\n+\nIIII\n")
        from gside.analysis.multigene_identifier import identify_multigene

        payload = identify_multigene(str(fq)).to_dict()
        assert payload["species"] == "Salmonella"
        assert payload["control_source"] == "fallback"
        assert payload["control_gene"] == "sample_median"
        assert payload["control_depth"] == 20.0
        assert any("sample median" in note.lower() for note in payload["notes"])


def _depth_hit(gene: str, depth: float, identity: float = 99.0, coverage: float = 100.0) -> dict:
    return {"gene": gene, "identity": identity, "coverage": coverage, "mean_depth": depth}


class TestControlDepth:
    def test_dedicated_control_uses_deepest_allele_per_gene(self):
        hits = [
            _depth_hit("lpxa", 3.0, identity=85.0),
            _depth_hit("lpxa", 30.0, identity=85.0),
        ]
        depth, label, source = mg._control_depth(hits, ["lpxa"])
        assert depth == pytest.approx(30.0)
        assert label == "lpxa"
        assert source == "dedicated"

    def test_control_below_identity_floor_is_ignored(self):
        depth, label, source = mg._control_depth(
            [_depth_hit("rpod", 30.0, identity=79.0)], ["rpod"]
        )
        assert depth is None
        assert label == ""
        assert source == ""

    def test_control_below_coverage_floor_is_ignored(self):
        depth, _label, _source = mg._control_depth(
            [_depth_hit("rpod", 30.0, coverage=49.0)], ["rpod"]
        )
        assert depth is None

    def test_fallback_is_sample_median_of_significant_markers(self):
        hits = [_depth_hit("inva", 10.0), _depth_hit("hilA", 30.0), _depth_hit("rpod", 50.0)]
        depth, label, source = mg._control_depth(hits, [])
        assert depth == pytest.approx(30.0)
        assert label == "sample_median"
        assert source == "fallback"

    def test_zero_depth_hits_are_ignored(self):
        depth, label, source = mg._control_depth([_depth_hit("rpod", 0.0)], ["rpod"])
        assert depth is None
        assert label == ""
        assert source == ""


class TestApplyDepthGate:
    def test_absolute_and_relative_floors(self):
        hits = [
            _depth_hit("rpod", 40.0),
            _depth_hit("inva", 30.0),
            _depth_hit("hilA", 2.0),
            _depth_hit("styM", 1.0),
        ]
        gate = mg._apply_depth_gate(hits, ["rpod"], 5.0, 0.10)
        assert [hit["gene"] for hit in gate.hits] == ["rpod", "inva"]
        assert {hit["gene"] for hit in gate.filtered} == {"hilA", "styM"}
        assert gate.control_gene == "rpod"
        assert gate.hits[1]["relative_depth_pct"] == pytest.approx(75.0)

    def test_missing_control_falls_back_to_absolute_floor_only(self):
        # No usable control at all (nothing passes the global identity floor):
        # the relative gate is disabled but the absolute floor still filters.
        hits = [
            _depth_hit("inva", 30.0, identity=80.0),
            _depth_hit("hilA", 2.0, identity=80.0),
        ]
        gate = mg._apply_depth_gate(hits, ["rpod"], 5.0, 0.10)
        assert [hit["gene"] for hit in gate.hits] == ["inva"]
        assert gate.control_depth is None

    def test_missing_dedicated_control_uses_sample_median(self):
        hits = [_depth_hit("inva", 30.0), _depth_hit("hilA", 2.0)]
        gate = mg._apply_depth_gate(hits, ["rpod"], 5.0, 0.10)
        assert gate.control_source == "fallback"
        assert gate.control_depth == 16.0
        assert [hit["gene"] for hit in gate.hits] == ["inva"]

    def test_assembly_hits_without_depth_pass_through(self):
        hits = [{"gene": "inva", "identity": 99.0, "coverage": 100.0}]
        gate = mg._apply_depth_gate(hits, ["rpod"], 5.0, 0.10)
        assert gate.hits == hits
        assert gate.control_depth is None

    def test_hit_without_depth_is_kept_ungated(self):
        hits = [
            {"gene": "inva", "identity": 99.0, "coverage": 100.0},
            _depth_hit("rpod", 40.0),
        ]
        gate = mg._apply_depth_gate(hits, ["rpod"], 5.0, 0.10)
        assert [hit["gene"] for hit in gate.hits] == ["inva", "rpod"]


class TestCustomMarkerDbInput:
    def _custom_db(self, tmp_path: Path) -> Path:
        root = tmp_path / "mydb"
        root.mkdir()
        (root / "markers.fasta").write_text(">markers~~~g1~~~x marker\nACGTACGT\n")
        (root / "marker_rules.yaml").write_text(
            "rules:\n- species: Test_species\n  genes: [g1]\n  min_hits: 1\n"
        )
        return root

    def test_string_marker_db_path_is_coerced(self, monkeypatch, tmp_path):
        from gside.analysis.multigene_identifier import identify_multigene
        from gside.marker_db import MarkerDatabase

        seen: dict[str, object] = {}

        def fake_blast(path, threads=4, marker_db=None):
            seen["db"] = marker_db
            return [{"gene": "g1", "identity": 99.0, "coverage": 100.0}]

        monkeypatch.setattr(mg, "_blast_contigs", fake_blast)
        genome = tmp_path / "g.fna"
        genome.write_text(">c1\nACGT\n")
        db_root = self._custom_db(tmp_path)
        result = identify_multigene(str(genome), marker_db=str(db_root))
        assert result.species == "Test_species"
        assert isinstance(seen["db"], MarkerDatabase)
        assert seen["db"].root == db_root


def test_main_text_output_prints_markers_and_notes(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr("sys.argv", ["multigene", str(tmp_path / "g.fna")])
    monkeypatch.setattr(
        mg,
        "_load_rules",
        lambda db=None: [
            {"species": "Shigella", "genes": ["uida", "ipah"], "min_hits": 2, "min_identity": 95}
        ],
    )
    monkeypatch.setattr(
        mg,
        "_blast_contigs",
        lambda path, threads=4, marker_db=None: [
            {"gene": "uida", "identity": 99.0, "coverage": 100.0},
            {"gene": "ipah", "identity": 99.5, "coverage": 100.0},
        ],
    )
    mg.main()
    out = capsys.readouterr().out
    assert "Species: Shigella" in out
    assert "uida: 99.0% identity" in out
    assert "NOTE:" in out


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
    def test_shigella_species_rules_upgrade_the_genus_call(self):
        # Genus signature alone stays at genus level; a species determinant
        # (rfc / wbgZ / sb1_wzx) upgrades the call via priority_over.
        assert _pick("uida", "ipah") == "Shigella"
        assert _pick("uida", "ipah", "rfc") == "Shigella_flexneri"
        assert _pick("uida", "ipah", "wbgz", "sonnei_psuk") == "Shigella_sonnei"
        assert _pick("uida", "ipah", "sb1_wzx") == "Shigella_boydii"
        # lacY-positive keeps the call E. coli: every Shigella rule is
        # lacY-negative (exclude_genes) so EIEC never fires Shigella.
        assert _pick("uida", "ipah", "lacy") == "Escherichia_coli"
        assert _pick("uida", "ipah", "rfc", "lacy") == "Escherichia_coli"

    def test_shigella_species_rules_ship(self):
        rules = {rule["species"]: rule for rule in _load_rules()}
        assert {
            "Shigella_flexneri",
            "Shigella_sonnei",
            "Shigella_dysenteriae",
            "Shigella_boydii",
        } <= set(rules)
        assert "Shigella" in rules
        for species in (
            "Shigella_flexneri",
            "Shigella_sonnei",
            "Shigella_dysenteriae",
            "Shigella_boydii",
        ):
            rule = rules[species]
            assert rule["priority_over"] == ["Shigella", "Escherichia_coli"]
            assert "lacy" in [gene.lower() for gene in rule.get("exclude_genes", [])]

    def test_shigella_call_notes_point_to_gapit(self, monkeypatch, tmp_path):
        monkeypatch.setattr(
            mg,
            "_blast_contigs",
            lambda path, threads: [
                {"gene": "uida", "identity": 99.0, "coverage": 100.0, "ref_len": 100},
                {"gene": "ipah", "identity": 100.0, "coverage": 100.0, "ref_len": 100},
            ],
        )
        result = mg.identify_multigene(str(tmp_path / "g.fna"))
        assert result.species == "Shigella"
        assert any("gapit" in note.lower() for note in result.notes)
        assert any("serotyp" in note.lower() for note in result.notes)

    def test_bordetella_pertussis_needs_both_targets(self):
        assert _pick("ptxs1", "tcfa") == "Bordetella_pertussis"
        assert _pick("ptxs1") == "Unknown"

    def test_renamed_labels(self):
        assert _pick("tlh") == "Vibrio_parahaemolyticus"
        assert _pick("inva") == "Salmonella_enterica"


def _bundled_marker_genes() -> set[str]:
    fasta = config.SPECIES_DB_DIR / "D1_marker" / "markers.fasta"
    return {
        line.split("~~~")[1].lower()
        for line in fasta.read_text(encoding="utf-8").splitlines()
        if line.startswith(">") and "~~~" in line
    }


class TestMarkerMetadata:
    def test_control_genes_exist_in_markers_fasta(self):
        controls, _defaults = mg._control_config()
        assert controls, "bundled db must ship control genes"
        assert set().union(*controls.values()) <= _bundled_marker_genes()

    def test_bundled_reads_defaults(self):
        _controls, defaults = mg._control_config()
        assert defaults == {"min_depth": 5.0, "min_relative_depth": 0.1, "min_mapq": 20}

    def test_missing_metadata_file_is_empty(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mg, "_MARKER_METADATA", tmp_path / "absent.yaml")
        assert mg._load_metadata() == {}
        assert mg._control_config() == ({}, {})

    def test_invalid_metadata_yaml_is_empty(self, monkeypatch, tmp_path):
        bad = tmp_path / "metadata.yaml"
        bad.write_text("control_genes: [unclosed", encoding="utf-8")
        monkeypatch.setattr(mg, "_MARKER_METADATA", bad)
        assert mg._load_metadata() == {}

    def test_non_dict_metadata_is_empty(self, monkeypatch, tmp_path):
        bad = tmp_path / "metadata.yaml"
        bad.write_text("- just\n- a\n- list\n", encoding="utf-8")
        monkeypatch.setattr(mg, "_MARKER_METADATA", bad)
        assert mg._load_metadata() == {}

    def test_control_config_accepts_string_gene(self, monkeypatch, tmp_path):
        meta = tmp_path / "metadata.yaml"
        meta.write_text("control_genes:\n  Salmonella: rpod\n", encoding="utf-8")
        monkeypatch.setattr(mg, "_MARKER_METADATA", meta)
        controls, _defaults = mg._control_config()
        assert controls == {"Salmonella": ["rpod"]}


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
        mg,
        "_blast_contigs",
        lambda path, threads: [
            {"gene": "b", "identity": 100.0, "coverage": 100.0, "ref_len": 100},
            {"gene": "a", "identity": 100.0, "coverage": 100.0, "ref_len": 100},
        ],
    )
    # No rule is required to exercise detected_markers ordering.
    result = mg.identify_multigene(str(tmp_path / "g.fna"))
    assert [m["gene"] for m in result.detected_markers] == ["a", "b"]
