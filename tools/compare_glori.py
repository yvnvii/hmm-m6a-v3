#!/usr/bin/env python
"""Compare hmm-m6a-v3 pileup counts with GLORI-tools' <prefix>.totalformat.txt.

For every (contig, 1-based position, strand), the v3 pileup's (n, A) is
compared with two sets of GLORI-tools counts:

  overall  the AGcov / Acov columns: every read, no A-count filter
           (matches a v3 pileup made with --max-unconverted -1)
  bin      the "3;total,AG,A" column: reads with <= 3 A's, what GLORI's
           caller tests at its default --cutoff 3 (matches v3
           --glori-compat / run-bam defaults). Other bins: --bin.

Each comparison reports the number of mismatching positions with examples.
GLORI-tools only reports positions covered by >= 2 reads, so a position
missing on its side counts as (0, 0). Contig suffixes such as _AG_converted
are stripped on both sides.

Optionally (--glori-bam with --v3-bam) the primary alignments of two BAMs
are compared as well, for a v3 run that did its own alignment.

Exits 0 when the comparison named by --expect (default: the bin) has no
mismatches (and the BAMs agree, if given), 1 otherwise.

usage:
  compare_glori.py --totalformat GLORI/x.totalformat.txt --pileup out/x.pileup.tsv
                   [--bin 3] [--expect bin|overall|none]
                   [--glori-bam GLORI/x_merged.sorted.bam --v3-bam out/x.bam]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter

import pandas as pd

SUFFIX = "_AG_converted"


def _strip(contig: str) -> str:
    return contig[: -len(SUFFIX)] if contig.endswith(SUFFIX) else contig


def bam_reads(path: str) -> Counter:
    import pysam

    c: Counter = Counter()
    with pysam.AlignmentFile(path) as fh:
        for a in fh:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            c[(_strip(a.reference_name), a.reference_start,
               "-" if a.is_reverse else "+", a.cigarstring, a.query_sequence)] += 1
    return c


def glori_counts(path: str, bin_: str = "3") -> tuple[dict, dict]:
    """(contig, 1-based pos, strand) -> (A+G, A), as (overall, bin) dicts."""
    overall, binned = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (_strip(f[0]), int(f[1]), f[2])
            overall[key] = (int(f[9]), int(f[10]))
            bins = dict(x.split(";", 1) for x in f[14:] if ";" in x)
            if bin_ not in bins:
                raise SystemExit(f"{path}: no '{bin_};' column at {key}; "
                                 f"bins present: {sorted(bins)}")
            _tot, ag, a = map(int, bins[bin_].split(","))
            binned[key] = (ag, a)
    return overall, binned


def v3_counts(path: str) -> dict:
    """(contig, 1-based pos, strand) -> (n, A) from a v3 pileup TSV."""
    df = pd.read_csv(path, sep="\t", dtype={"contig": str})
    out = {}
    for r in df.itertuples():
        pos0 = (r.region_start + r.pos0_in_region if r.strand == "+"
                else r.region_end - 1 - r.pos0_in_region)
        out[(_strip(r.contig), pos0 + 1, r.strand)] = (int(r.n), int(r.A))
    return out


def diff_counts(g: dict, v: dict) -> tuple[list, list]:
    """Positions v3 scored whose counts differ from GLORI's, and GLORI
    positions with coverage that v3 did not score."""
    mismatch = [(k, g.get(k, (0, 0)), v[k]) for k in sorted(v)
                if g.get(k, (0, 0)) != v[k]]
    missing = sorted(k for k in g if k not in v and g[k][0] > 0)
    return mismatch, missing


def _report(label: str, g: dict, v: dict, show: int) -> bool:
    mismatch, missing = diff_counts(g, v)
    print(f"{label}: positions GLORI {len(g)}  v3 {len(v)}  "
          f"mismatching {len(mismatch)}  GLORI-only with coverage {len(missing)}")
    for (c, p, s), gc, vc in mismatch[:show]:
        print(f"    {c}:{p}({s})  GLORI(AG,A)={gc}  v3(n,A)={vc}")
    for c, p, s in missing[:show]:
        print(f"    {c}:{p}({s})  GLORI(AG,A)={g[(c, p, s)]}  not in v3 pileup")
    return not mismatch and not missing


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("\n\n", 1)[1])
    p.add_argument("--totalformat", "--glori-format", required=True,
                   help="GLORI-tools <prefix>.totalformat.txt")
    p.add_argument("--pileup", "--v3-pileup", required=True, help="hmm-m6a-v3 pileup TSV")
    p.add_argument("--bin", default="3",
                   help="GLORI A-count bin to compare (default 3, its caller's --cutoff)")
    p.add_argument("--expect", choices=["bin", "overall", "none"], default="bin",
                   help="comparison that must match for exit 0 (default bin; "
                        "overall for a --max-unconverted -1 pileup; none = report only)")
    p.add_argument("--glori-bam", help="GLORI-tools <prefix>_merged.sorted.bam (optional)")
    p.add_argument("--v3-bam", help="hmm-m6a-v3 BAM (optional, with --glori-bam)")
    p.add_argument("--show", type=int, default=10, help="examples to print per comparison")
    a = p.parse_args(argv)
    if bool(a.glori_bam) != bool(a.v3_bam):
        p.error("--glori-bam and --v3-bam go together")

    ok = True
    if a.glori_bam:
        rg, rv = bam_reads(a.glori_bam), bam_reads(a.v3_bam)
        only_g, only_v = rg - rv, rv - rg
        print(f"reads: GLORI {sum(rg.values())}  v3 {sum(rv.values())}  "
              f"shared {sum((rg & rv).values())}  only-GLORI {sum(only_g.values())}  "
              f"only-v3 {sum(only_v.values())}")
        for k in list(only_g)[:a.show]:
            print("    only-GLORI", k[:4], k[4][:40])
        for k in list(only_v)[:a.show]:
            print("    only-v3   ", k[:4], k[4][:40])
        ok &= not only_g and not only_v

    overall, binned = glori_counts(a.totalformat, a.bin)
    v = v3_counts(a.pileup)
    ok_overall = _report("overall (AGcov, Acov)", overall, v, a.show)
    ok_bin = _report(f"bin {a.bin} (AG, A)     ", binned, v, a.show)
    if a.expect == "bin":
        ok &= ok_bin
    elif a.expect == "overall":
        ok &= ok_overall
    print(f"expected match ({a.expect}): " + ("MATCH" if ok else "DIFFER")
          if a.expect != "none" else "report only")
    return 0 if ok or a.expect == "none" else 1


if __name__ == "__main__":
    sys.exit(main())
