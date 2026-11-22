"""Tests for the FASTQ -> calls pipeline on GLORI-tools style references.

Unit tests run anywhere. The end-to-end tests need STAR and bowtie on PATH
(they build a tiny synthetic reference in a temp dir) and are skipped
otherwise.
"""
from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

from hmm_m6a_fastq.align import (AlignStats, GloriRefs, _process_stage, encode_name,
                                  iter_fastq, restore_bases, restore_bases_mirrored)
from hmm_m6a_fastq.annotation import TranscriptModel, annotate_sites, load_tbl2
from hmm_m6a_fastq.cli import main
from hmm_m6a_fastq.regions import _merge

EXAMPLE = Path(__file__).resolve().parents[1] / "example"
have_tools = shutil.which("STAR") and shutil.which("bowtie") and shutil.which("bowtie-build")
needs_tools = pytest.mark.skipif(not have_tools, reason="STAR / bowtie not installed")


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("seq", ["ACGTAAN", "GGG", "A" * 150, "TTTTA"])
def test_name_mask_round_trip(seq):
    name = encode_name(7, seq)
    assert name.startswith("7_")
    assert restore_bases(seq.replace("A", "G"), name) == seq


def test_mirrored_restore_matches_glori_tools():
    """GLORI-tools restores rvsCom forward hits via revcomp -> put A -> revcomp."""
    rc = str.maketrans("ACGT", "TGCA")
    seq = "ACGTAAGGCA"
    name = encode_name(0, seq)
    conv = seq.replace("A", "G")
    glori = list(conv.translate(rc)[::-1])
    for k, b in enumerate(seq):
        if b == "A":
            glori[k] = "A"
    glori = "".join(glori).translate(rc)[::-1]
    assert restore_bases_mirrored(conv, name) == glori


def _stage_sam(tmp_path: Path, recs):
    """recs: list of (name, seq, flag). Writes a SAM on one 100-nt contig."""
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c", "LN": 100}]}
    p = tmp_path / "stage.sam"
    with pysam.AlignmentFile(str(p), "w", header=hdr) as out:
        for name, seq, flag in recs:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.query_sequence = name, flag, seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            if not flag & 4:
                a.reference_id, a.reference_start = 0, 0
                a.mapping_quality, a.cigartuples = 255, [(0, len(seq))]
            out.write(a)
    return p


@pytest.mark.parametrize("stage", ["ag", "rvs"])
@pytest.mark.parametrize("compat", [False, True])
def test_stage_strand_rules(tmp_path: Path, stage, compat):
    """Default: keep the stage's own strand, retry everything else.
    glori_compat (GLORI-tools): stage 1 drops reverse hits instead of retrying
    them, stage 2 keeps both orientations; only unmapped reads go on."""
    seq = "ACGTAAGGCATTGCA"
    fwd, rev, un = (encode_name(i, seq) for i in range(3))
    conv = seq.replace("A", "G")
    rc = str.maketrans("ACGT", "TGCA")
    revcomp = lambda x: x.translate(rc)[::-1]          # SAM stores reverse hits revcomp'd
    sam = _stage_sam(tmp_path, [(fwd, conv, 0), (rev, revcomp(conv), 16), (un, conv, 4)])
    hdr = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "c", "LN": 100}]})
    stats = AlignStats()
    nxt = tmp_path / "next.fq"
    with pysam.AlignmentFile(str(tmp_path / "o.bam"), "wb", header=hdr) as bam, \
            open(nxt, "w") as nf:
        _process_stage(sam, stage, bam, nf, {"c": 0}, stats, glori_compat=compat)
    passed = [l[1:].strip() for l in nxt.read_text().splitlines()[::4]]
    kept = {a.query_name: a for a in pysam.AlignmentFile(str(tmp_path / "o.bam"))}
    if not compat:
        own = fwd if stage == "ag" else rev
        assert set(kept) == {own}
        assert kept[own].query_sequence == (seq if own == fwd else revcomp(seq))
        assert sorted(passed) == sorted({fwd, rev, un} - {own})
        assert stats.wrong_strand == 0
    elif stage == "ag":
        assert passed == [un]
        assert set(kept) == {fwd} and stats.wrong_strand == 1
        assert kept[fwd].query_sequence == seq
    else:
        assert passed == [un]
        assert set(kept) == {fwd, rev} and stats.wrong_strand == 0
        assert kept[fwd].query_sequence == restore_bases_mirrored(conv, fwd)
        assert kept[rev].query_sequence == revcomp(seq)


@pytest.mark.parametrize("argv, expect", [
    ([], {"min_mapq": 20, "min_baseq": 20, "max_unconverted": -1}),
    (["--glori-compat"], {"min_mapq": 0, "min_baseq": 10, "max_unconverted": 3}),
    (["--glori-compat", "--max-unconverted", "-1"],
     {"min_mapq": 0, "min_baseq": 10, "max_unconverted": -1}),
    (["--max-unconverted", "3"], {"min_mapq": 20, "min_baseq": 20, "max_unconverted": 3}),
])
def test_pileup_filter_defaults(argv, expect):
    """--glori-compat switches the filter defaults; explicit values win."""
    from hmm_m6a_fastq.cli import _pileup_filters, build_parser
    args = build_parser().parse_args(["pileup", "--bam", "x.bam", "--genome", "g.fa",
                                      "-o", "o.tsv", *argv])
    assert _pileup_filters(args) == expect


def test_iter_fastq_plain_and_malformed(tmp_path: Path):
    p = tmp_path / "x.fq"
    p.write_text("@r1 extra\nacgt\n+\nIIII\n@r2\nGG\n+\nII\n")
    assert list(iter_fastq([str(p)])) == [("r1", "ACGT", "IIII"), ("r2", "GG", "II")]
    p.write_text("@r1\nACGT\n+\nII\n")
    with pytest.raises(ValueError):
        list(iter_fastq([str(p)]))


def test_lift_plus_strand_spliced():
    m = TranscriptModel("T", "chr1", "+", "G", [(100, 110), (200, 220)])
    assert m.lift(0, 5) == (100, [(0, 5)])
    assert m.lift(8, 6) == (108, [(0, 2), (3, 90), (0, 4)])
    assert m.lift(25, 6) is None                      # runs off the 30-nt transcript


def test_lift_minus_strand_spliced():
    # transcript 5' end is genomic 219; exon 2 (200-220) comes first
    m = TranscriptModel("T", "chr1", "-", "G", [(100, 110), (200, 220)])
    assert m.lift(0, 5) == (215, [(0, 5)])
    assert m.lift(18, 6) == (106, [(0, 4), (3, 90), (0, 2)])


def test_load_tbl2(tmp_path: Path):
    p = tmp_path / "a.tbl2"
    p.write_text("0\tNM_1\tchr2\t-\t5\t50\t5\t50\t2\t5,30,\t10,50,\t0\tGENEX\tcmpl\tcmpl\t0,\n")
    m = load_tbl2(str(p))["NM_1"]
    assert (m.chrom, m.strand, m.gene, m.exons, m.length) == ("chr2", "-", "GENEX", [(5, 10), (30, 50)], 25)


def test_merge_clusters():
    s, e, n = _merge(np.array([0, 10, 100, 5]), np.array([20, 30, 120, 8]), max_gap=50)
    assert s.tolist() == [0, 100] and e.tolist() == [30, 120] and n.tolist() == [3, 1]
    s, e, n = _merge(np.array([0, 100]), np.array([20, 120]), max_gap=80)
    assert s.tolist() == [0] and e.tolist() == [120]


def test_annotate_sites(tmp_path: Path):
    base = tmp_path / "b.base"
    base.write_text("chr1\t10\t11\t+\tG1\tNM_1\tGeneOne\tNM_1\tprotein_coding\n"
                    "chr1\t10\t11\t-\tG2\tNM_2\tGeneTwo\tNM_2\tlncRNA\n")
    df = pd.DataFrame({"contig": ["chr1", "chr1", "chr2"], "pos0": [10, 10, 10],
                       "strand": ["+", "-", "+"], "call": [True] * 3})
    out = annotate_sites(df, str(base), quiet=True)
    assert out["gene"].tolist() == ["G1", "G2", "ELSE"]
    assert out["biotype"].tolist() == ["protein_coding", "lncRNA", "NA"]
    assert annotate_sites(df, None)["gene"].tolist() == ["ELSE"] * 3


def test_refs_report_missing_files(tmp_path: Path):
    (tmp_path / "index/genome").mkdir(parents=True)
    with pytest.raises(ValueError):
        GloriRefs.resolve(str(tmp_path))


# ---------------------------------------------------------------------------
# End-to-end on the synthetic GLORI-tools layout
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def dataset(tmp_path_factory) -> Path:
    if not have_tools:
        pytest.skip("STAR / bowtie not installed")
    spec = importlib.util.spec_from_file_location("build_example", EXAMPLE / "build_example.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    d = tmp_path_factory.mktemp("glori")
    mod.build_indices(mod.write_dataset(d), threads=1)
    return d


def _truth(d: Path) -> set:
    t = pd.read_csv(d / "truth.bed", sep="\t", header=None)
    return set(zip(t[0], t[1]))


def _called(tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t")
    return df[df["call"]]


@needs_tools
def test_refs_resolve_layout(dataset: Path):
    r = GloriRefs.resolve(str(dataset / "ref"))
    assert r.genome.endswith("syn.main.fa")
    assert GloriRefs.star_dir(r.genome_ag).endswith("index/genome/syn.AG_conversion")
    assert r.baseanno.endswith(".noredundance.base")


@needs_tools
def test_run_end_to_end(dataset: Path, tmp_path: Path):
    rc = main(["-q", "run", "--fastq", str(dataset / "full.clean.fq"),
               "--ref-dir", str(dataset / "ref"), "--outdir", str(tmp_path), "--threads", "2"])
    assert rc == 0
    stats = json.loads((tmp_path / "sample.bam.stats.json").read_text())
    assert stats["unaligned"] == 50                      # exactly the junk reads
    assert stats["genome_plus"] > 0 and stats["genome_minus"] > 0
    assert stats["transcriptome"] + stats["genome_plus"] + stats["genome_minus"] == stats["reads_in"] - 50

    called = _called(tmp_path / "sample.calls.tsv")
    truth = _truth(dataset)
    got = set(zip(called["contig"], called["pos0"]))
    assert len(truth & got) >= 0.9 * len(truth)
    assert len(got - truth) == 0
    assert (called["gene"] != "ELSE").all()


@needs_tools
def test_bam_bases_match_genome(dataset: Path, tmp_path: Path):
    bam = tmp_path / "s.bam"
    assert main(["-q", "align", "--fastq", str(dataset / "full.clean.fq"),
                 "--ref-dir", str(dataset / "ref"), "-o", str(bam), "--threads", "2"]) == 0
    fa = pysam.FastaFile(str(dataset / "ref/index/genome/syn.main.fa"))
    checked = bad = 0
    with pysam.AlignmentFile(str(bam)) as fh:
        for r in fh:
            ref = fa.fetch(r.reference_name, r.reference_start, r.reference_end)
            tgt, conv = ("T", "C") if r.is_reverse else ("A", "G")
            for qp, rp in r.get_aligned_pairs(matches_only=True):
                rb, qb = ref[rp - r.reference_start], r.query_sequence[qp]
                checked += 1
                ok = qb in (tgt, conv) if rb == tgt else qb == rb
                bad += not ok
    # error-free simulation: only STAR's junction-edge placements may disagree
    assert checked > 100_000 and bad <= 10


@needs_tools
def test_two_stage_train_then_call(dataset: Path, tmp_path: Path):
    ref = str(dataset / "ref")
    bam, pile = tmp_path / "s.bam", tmp_path / "s.pileup.tsv"
    params, calls = tmp_path / "p.json", tmp_path / "c.tsv"
    assert main(["-q", "align", "--fastq", str(dataset / "full.clean.fq"),
                 "--ref-dir", ref, "-o", str(bam)]) == 0
    assert main(["-q", "pileup", "--bam", str(bam), "--ref-dir", ref, "-o", str(pile)]) == 0
    assert main(["-q", "train", "--pileup", str(pile), "-o", str(params)]) == 0
    assert main(["-q", "call", "--pileup", str(pile), "--params", str(params),
                 "--ref-dir", ref, "-o", str(calls), "--threads", "1"]) == 0
    got = _called(calls)
    assert len(_truth(dataset) & set(zip(got["contig"], got["pos0"]))) >= 30
