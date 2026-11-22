"""scripts/compare_glori.py: parsing and diffing of GLORI totalformat vs v3 pileup."""
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "compare_glori.py"
spec = importlib.util.spec_from_file_location("compare_glori", SCRIPT)
cg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cg)


def _glori_line(contig, pos1, strand, ag_all, a_all, bin3):
    f = [contig, str(pos1), strand] + ["x"] * 6 + [str(ag_all), str(a_all)] + ["x"] * 3
    return "\t".join(f + [f"3;{bin3[0]},{bin3[1]},{bin3[2]}", "all;0,0,0"]) + "\n"


def _write(tmp_path):
    g = tmp_path / "g.totalformat.txt"
    g.write_text(_glori_line("chr1", 11, "+", 10, 4, (9, 8, 3))
                 + _glori_line("chr1", 20, "-", 5, 1, (5, 5, 1)))
    v = tmp_path / "v.tsv"
    # + strand: pos0 = region_start + pos0_in_region; - strand: region_end - 1 - pos0_in_region
    v.write_text("contig\tstrand\tregion_start\tregion_end\tpos0_in_region\tn\tA\n"
                 "chr1\t+\t5\t30\t5\t8\t3\n"
                 "chr1\t-\t5\t30\t10\t5\t1\n"
                 "chr1\t+\t5\t30\t20\t1\t0\n")
    return g, v


def test_counts_match_at_cutoff_bin(tmp_path):
    g, v = _write(tmp_path)
    gc, vc = cg.glori_counts(str(g), "3"), cg.v3_counts(str(v))
    assert gc[("chr1", 11, "+")] == (8, 3) and vc[("chr1", 11, "+")] == (8, 3)
    assert vc[("chr1", 20, "-")] == (5, 1)
    mismatch, missing = cg.diff_counts(gc, vc)
    # chr1:26 + has 1 read in v3 but GLORI drops it (< 2 reads) -> (0, 0) mismatch
    assert [m[0] for m in mismatch] == [("chr1", 26, "+")] and missing == []


def test_all_bin_and_exit_code(tmp_path):
    g, v = _write(tmp_path)
    gc = cg.glori_counts(str(g), "all")
    assert gc[("chr1", 11, "+")] == (10, 4)
    assert cg.main(["--glori-format", str(g), "--v3-pileup", str(v)]) == 1
