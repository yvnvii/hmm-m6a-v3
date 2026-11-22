#!/usr/bin/env python
"""Check that hmm-m6a-v3 and GLORI-tools used the same reads and counts.

Compares a v3 run (normally made with --glori-compat) against GLORI-tools
on the same FASTQ and reference:

1. Reads: the multiset of (contig, start, strand, CIGAR, sequence) of the
   primary alignments in each BAM. Skipped unless both BAMs are given.
2. Counts: GLORI's caller input (A+G, A) per (contig, 1-based pos, strand)
   from its *totalformat.txt, against v3's (n, A) from the pileup TSV.
   By default GLORI's bin for reads with <= 3 A's is used (its caller's
   --cutoff 3), matching v3's --max-unconverted 3; pass --acut all for its
   unfiltered AGcov/Acov columns (v3 --max-unconverted -1). GLORI only
   reports positions covered by >= 2 reads, so a position missing on its
   side counts as (0, 0).

Exits 0 when reads (if compared) and counts all agree, 1 otherwise.

usage:
  compare_glori.py --glori-format GLORI.totalformat.txt --v3-pileup v3.pileup.tsv
                   [--glori-bam GLORI.bam --v3-bam v3.bam] [--acut 3|all]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter

import pandas as pd


def bam_reads(path: str) -> Counter:
    import pysam

    c: Counter = Counter()
    with pysam.AlignmentFile(path) as fh:
        for a in fh:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            c[(a.reference_name.replace("_AG_converted", ""), a.reference_start,
               "-" if a.is_reverse else "+", a.cigarstring, a.query_sequence)] += 1
    return c


def glori_counts(path: str, acut: str = "3") -> dict:
    """(contig, 1-based pos, strand) -> (A+G, A) from GLORI's totalformat."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if acut == "all":
                ag, a = int(f[9]), int(f[10])
            else:
                bins = dict(x.split(";", 1) for x in f[14:])
                _tot, ag, a = map(int, bins[acut].split(","))
            out[(f[0], int(f[1]), f[2])] = (ag, a)
    return out


def v3_counts(path: str) -> dict:
    """(contig, 1-based pos, strand) -> (n, A) from a v3 pileup TSV."""
    df = pd.read_csv(path, sep="\t")
    out = {}
    for r in df.itertuples():
        pos0 = (r.region_start + r.pos0_in_region if r.strand == "+"
                else r.region_end - 1 - r.pos0_in_region)
        out[(r.contig, pos0 + 1, r.strand)] = (int(r.n), int(r.A))
    return out


def diff_counts(g: dict, v: dict) -> tuple[list, list]:
    """Positions v3 scored whose counts differ from GLORI's, and GLORI
    positions with coverage that v3 did not score."""
    mismatch = [(k, g.get(k, (0, 0)), v[k]) for k in sorted(v)
                if g.get(k, (0, 0)) != v[k]]
    missing = sorted(k for k in g if k not in v and g[k][0] > 0)
    return mismatch, missing


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--glori-format", required=True, help="GLORI-tools *totalformat.txt")
    p.add_argument("--v3-pileup", required=True, help="hmm-m6a-v3 pileup TSV")
    p.add_argument("--glori-bam", help="GLORI-tools final BAM (optional)")
    p.add_argument("--v3-bam", help="hmm-m6a-v3 BAM (optional)")
    p.add_argument("--acut", default="3",
                   help="GLORI A-count bin to compare, or 'all' for no read filter (default 3)")
    p.add_argument("--show", type=int, default=10, help="examples to print per category")
    a = p.parse_args(argv)

    ok = True
    if bool(a.glori_bam) != bool(a.v3_bam):
        p.error("--glori-bam and --v3-bam go together")
    if a.glori_bam:
        rg, rv = bam_reads(a.glori_bam), bam_reads(a.v3_bam)
        only_g, only_v = rg - rv, rv - rg
        print(f"reads: GLORI {sum(rg.values())}  v3 {sum(rv.values())}  "
              f"shared {sum((rg & rv).values())}  only-GLORI {sum(only_g.values())}  "
              f"only-v3 {sum(only_v.values())}")
        for k in list(only_g)[:a.show]:
            print("  only-GLORI", k[:4], k[4][:40])
        for k in list(only_v)[:a.show]:
            print("  only-v3   ", k[:4], k[4][:40])
        ok &= not only_g and not only_v

    g, v = glori_counts(a.glori_format, a.acut), v3_counts(a.v3_pileup)
    mismatch, missing = diff_counts(g, v)
    print(f"sites: GLORI {len(g)}  v3 {len(v)}  count mismatches {len(mismatch)}  "
          f"GLORI sites with n>0 not scored by v3 {len(missing)}")
    for k, gc, vc in mismatch[:a.show]:
        print(f"  {k}  GLORI(AG,A)={gc}  v3(n,A)={vc}")
    for k in missing[:a.show]:
        print(f"  not in v3: {k}  GLORI(AG,A)={g[k]}")
    ok &= not mismatch and not missing
    print("MATCH" if ok else "DIFFER")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
