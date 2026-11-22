"""tools/compare_glori.py: GLORI totalformat vs v3 pileup counts."""
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "compare_glori.py"
spec = importlib.util.spec_from_file_location("compare_glori", SCRIPT)
cg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cg)


def _glori_line(contig, pos1, strand, overall, bin3):
    f = [contig, str(pos1), strand] + ["x"] * 6 + [str(overall[0]), str(overall[1])] + ["x"] * 3
    return "\t".join(f + ["2;0,0,0", "3;{},{},{}".format(*bin3)]) + "\n"


def _write(tmp_path, v_rows):
    g = tmp_path / "g.totalformat.txt"
    g.write_text(_glori_line("chr1_AG_converted", 11, "+", (10, 4), (9, 8, 3))
                 + _glori_line("chr1", 20, "-", (5, 1), (5, 5, 1)))
    v = tmp_path / "v.tsv"
    # + strand: pos0 = region_start + pos0_in_region; - strand: region_end - 1 - pos0_in_region
    v.write_text("contig\tstrand\tregion_start\tregion_end\tpos0_in_region\tn\tA\n" + v_rows)
    return g, v


FILTERED = "chr1\t+\t5\t30\t5\t8\t3\nchr1\t-\t5\t30\t10\t5\t1\n"
UNFILTERED = "chr1\t+\t5\t30\t5\t10\t4\nchr1\t-\t5\t30\t10\t5\t1\n"


def test_parse_overall_and_bin(tmp_path):
    g, v = _write(tmp_path, FILTERED)
    overall, binned = cg.glori_counts(str(g), "3")
    assert overall[("chr1", 11, "+")] == (10, 4)          # suffix stripped
    assert binned[("chr1", 11, "+")] == (8, 3)
    assert cg.v3_counts(str(v)) == {("chr1", 11, "+"): (8, 3), ("chr1", 20, "-"): (5, 1)}


def test_expected_comparison_decides_exit(tmp_path, capsys):
    g, v = _write(tmp_path, FILTERED)
    base = ["--totalformat", str(g), "--pileup", str(v)]
    assert cg.main(base) == 0
    out = capsys.readouterr().out
    assert "overall (AGcov, Acov): positions GLORI 2  v3 2  mismatching 1" in out
    assert "chr1:11(+)  GLORI(AG,A)=(10, 4)  v3(n,A)=(8, 3)" in out
    assert cg.main(base + ["--expect", "overall"]) == 1
    g, v = _write(tmp_path, UNFILTERED)
    assert cg.main(base + ["--expect", "overall"]) == 0
    assert cg.main(base) == 1


def test_single_read_positions_count_as_zero(tmp_path):
    g, v = _write(tmp_path, FILTERED + "chr1\t+\t5\t30\t20\t1\t0\n")
    _, binned = cg.glori_counts(str(g), "3")
    mismatch, missing = cg.diff_counts(binned, cg.v3_counts(str(v)))
    assert [m[0] for m in mismatch] == [("chr1", 26, "+")] and missing == []
