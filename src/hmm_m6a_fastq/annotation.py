"""GLORI-tools annotation files.

tbl2  (`..._genomic.gtf_change2Ens.tbl2`)
    UCSC genePred-with-bin table, one transcript per line:
        0 bin  1 transcript  2 chrom  3 strand  4 txStart  5 txEnd
        6 cdsStart  7 cdsEnd  8 exonCount  9 exonStarts  10 exonEnds
        11 score  12 gene ...
    exonStarts / exonEnds are comma-terminated, 0-based half-open.
    Used to lift transcriptome alignments back to genome coordinates.

noredundance.base  (`....tbl2.noredundance.base`)
    One line per exonic base, at most one gene per (chrom, pos, strand):
        0 chrom  1 pos0  2 pos1  3 strand  4 gene_id  5 transcript
        6 gene_name  7 isoform_name  8 biotype
    Used to label each site with its gene; unannotated sites get "ELSE".
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

TX_SUFFIX = "_AG_converted"


@dataclass
class TranscriptModel:
    transcript: str
    chrom: str
    strand: str
    gene: str
    exons: List[Tuple[int, int]]   # genomic, 0-based half-open, ascending

    @property
    def length(self) -> int:
        return sum(e - s for s, e in self.exons)

    def lift(self, tpos: int, length: int) -> Optional[Tuple[int, List[Tuple[int, int]]]]:
        """Map transcript interval [tpos, tpos+length) to the genome.

        `tpos` is on the transcript (5'->3'). Returns (genomic_start, cigar)
        with M blocks over exons and N blocks over introns, ascending genome
        order; None if the interval runs off the transcript.
        """
        L = self.length
        if tpos < 0 or length <= 0 or tpos + length > L:
            return None
        # Position in the exon concatenation taken in ascending genome order.
        a = tpos if self.strand == "+" else L - tpos - length
        b = a + length
        segs: List[Tuple[int, int]] = []
        off = 0
        for s, e in self.exons:
            n = e - s
            lo, hi = max(a, off), min(b, off + n)
            if lo < hi:
                segs.append((s + lo - off, s + hi - off))
            off += n
            if off >= b:
                break
        cigar: List[Tuple[int, int]] = []
        for k, (s, e) in enumerate(segs):
            if k:
                cigar.append((3, s - segs[k - 1][1]))   # N
            cigar.append((0, e - s))                      # M
        return segs[0][0], cigar


def load_tbl2(path: str) -> Dict[str, TranscriptModel]:
    """Read a GLORI-tools tbl2 file into {transcript_id: TranscriptModel}.

    If a transcript id occurs more than once (e.g. PAR genes on chrX/chrY),
    the last line wins, as in GLORI-tools' Transcrip2genome.py.
    """
    out: Dict[str, TranscriptModel] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            starts = [int(x) for x in f[9].split(",") if x]
            ends = [int(x) for x in f[10].split(",") if x]
            exons = sorted(zip(starts, ends))
            out[f[1]] = TranscriptModel(
                transcript=f[1], chrom=f[2], strand=f[3], gene=f[12], exons=exons)
    return out


BASE_COLS = ["gene", "gene_name", "transcript", "biotype"]


def annotate_sites(
    df: pd.DataFrame,
    base_path: Optional[str],
    chunksize: int = 5_000_000,
    quiet: bool = False,
) -> pd.DataFrame:
    """Add gene / gene_name / transcript / biotype columns to a calls table.

    Matches on (contig, strand, pos0). Streams the base file in chunks, so
    the genome-wide file never has to fit in memory. Sites with no entry get
    gene = "ELSE" (GLORI-tools convention); without a base file, all do.
    """
    df = df.copy()
    df["gene"] = "ELSE"
    for c in BASE_COLS[1:]:
        df[c] = "NA"
    if base_path is None or df.empty:
        return df

    # Query positions per (chrom, strand), as sorted int arrays.
    want: Dict[Tuple[str, str], np.ndarray] = {
        k: np.unique(g["pos0"].to_numpy(np.int64))
        for k, g in df.groupby(["contig", "strand"], sort=False)
    }
    hits = []
    reader = pd.read_csv(
        base_path, sep="\t", header=None, usecols=[0, 1, 3, 4, 5, 6, 8],
        names=["contig", "pos0", "strand", "gene", "transcript", "gene_name", "biotype"],
        dtype={0: str, 1: np.int64, 3: str, 4: str, 5: str, 6: str, 8: str},
        chunksize=chunksize,
    )
    n_lines = 0
    for chunk in reader:
        n_lines += len(chunk)
        for key, g in chunk.groupby(["contig", "strand"], sort=False):
            q = want.get(key)
            if q is None:
                continue
            m = np.isin(g["pos0"].to_numpy(np.int64), q, assume_unique=False)
            if m.any():
                hits.append(g[m])
    if not quiet:
        print(f"[annotate] scanned {n_lines} annotated bases", file=sys.stderr)
    if not hits:
        return df
    ann = pd.concat(hits).drop_duplicates(["contig", "strand", "pos0"], keep="first")
    merged = df.drop(columns=BASE_COLS).merge(ann, on=["contig", "strand", "pos0"], how="left")
    merged["gene"] = merged["gene"].fillna("ELSE")
    for c in BASE_COLS[1:]:
        merged[c] = merged[c].fillna("NA")
    return merged
