"""Where to run HMM chains, and running pileup / calling in parallel.

Without a regions BED, chains are built from the reads themselves: on each
strand, aligned blocks are merged into clusters whenever they are within
`max_gap` nt of each other, and clusters with fewer than `min_reads` blocks
are dropped (they cannot reach the calling depth anyway). Spliced reads
contribute one block per exon, so introns split clusters. Positions between
clusters carry no data, and the per-nucleotide HMM would treat them as
neutral emissions, so splitting there changes almost nothing while keeping
the work proportional to the covered territory instead of the genome.
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pysam

from hmm_m6a_fastq.calling import call_sites
from hmm_m6a_fastq.hmm import HMMParams
from hmm_m6a_fastq.pileup import Region, pileup_bam
from hmm_m6a_fastq.train import TranscriptObs

LIBRARY = "single-forward"   # GLORI reads are the RNA sense strand
# GLORI-tools' merged BAM keeps the converted-genome contig names.
GLORI_CONTIG_SUFFIX = "_AG_converted"


def strip_contig(name: str, suffix: Optional[str]) -> str:
    """`name` without a trailing `suffix` (the reference FASTA / output name)."""
    if suffix and name.endswith(suffix):
        return name[: -len(suffix)]
    return name


def with_bam_contigs(
    regions: Sequence[Region], bam_path: str, suffix: Optional[str]
) -> List[Region]:
    """Point regions named by reference contig (e.g. from a BED) at the BAM's
    contig names, for a BAM whose names carry `suffix`."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        to_bam = {strip_contig(c, suffix): c for c in bam.references}
    return [Region(name=r.name, contig=r.contig, start=r.start, end=r.end,
                   strand=r.strand,
                   bam_contig=None if to_bam.get(r.contig, r.contig) == r.contig
                   else to_bam[r.contig])
            for r in regions]


def _merge(starts: np.ndarray, ends: np.ndarray, max_gap: int):
    order = np.argsort(starts, kind="stable")
    s, e = starts[order], ends[order]
    run_end = np.maximum.accumulate(e)
    new = np.empty(len(s), dtype=bool)
    new[0] = True
    new[1:] = s[1:] > run_end[:-1] + max_gap
    cid = np.cumsum(new) - 1
    first = np.flatnonzero(new)
    last = np.append(first[1:] - 1, len(s) - 1)
    return s[first], run_end[last], np.bincount(cid)


def regions_from_coverage(
    bam_path: str,
    min_reads: int = 10,
    max_gap: int = 50,
    min_mapq: int = 20,
    strip_suffix: Optional[str] = None,
) -> List[Region]:
    """Strand-specific read clusters, one Region (HMM chain) each.

    With `strip_suffix`, region contigs (and names) drop that suffix from the
    BAM's contig names; `Region.bam_contig` keeps the BAM name for reading.
    """
    out: List[Region] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for bam_contig in bam.references:
            contig = strip_contig(bam_contig, strip_suffix)
            blocks = {"+": ([], []), "-": ([], [])}
            for aln in bam.fetch(bam_contig):
                if (aln.is_unmapped or aln.is_secondary or aln.is_supplementary
                        or aln.mapping_quality < min_mapq):
                    continue
                strand = "-" if aln.is_reverse else "+"
                bs, be = blocks[strand]
                for s, e in aln.get_blocks():
                    bs.append(s)
                    be.append(e)
            for strand, (bs, be) in blocks.items():
                if not bs:
                    continue
                cs, ce, n = _merge(np.asarray(bs), np.asarray(be), max_gap)
                for s, e, k in zip(cs, ce, n):
                    if k >= min_reads:
                        out.append(Region(
                            name=f"{contig}:{s + 1}-{e}({strand})",
                            contig=contig, start=int(s), end=int(e), strand=strand,
                            bam_contig=None if bam_contig == contig else bam_contig))
    return out


def regions_to_bed(regions: Sequence[Region], path: str) -> None:
    with open(path, "w") as fh:
        for r in regions:
            fh.write(f"{r.contig}\t{r.start}\t{r.end}\t{r.name}\t0\t{r.strand}\n")


# ---------------------------------------------------------------------------
# Parallel helpers (worker functions must be top-level to be picklable)
# ---------------------------------------------------------------------------
def _chunks(items: Sequence, n: int) -> List[Sequence]:
    if not items:
        return []
    size = max(1, -(-len(items) // n))
    return [items[i:i + size] for i in range(0, len(items), size)]


def _pileup_worker(args):
    bam, fasta, regions, kw = args
    return pileup_bam(bam, fasta, regions, show_progress=False, **kw)


def parallel_pileup(
    bam: str, fasta: str, regions: Sequence[Region], threads: int = 1, **kw
) -> List[Tuple[Region, TranscriptObs]]:
    kw.setdefault("library", LIBRARY)
    if threads <= 1 or len(regions) < 2:
        return pileup_bam(bam, fasta, regions, show_progress=False, **kw)
    jobs = [(bam, fasta, c, kw) for c in _chunks(list(regions), threads * 4)]
    out: List[Tuple[Region, TranscriptObs]] = []
    with ProcessPoolExecutor(max_workers=threads) as ex:
        for res in ex.map(_pileup_worker, jobs):
            out.extend(res)
    return out


def _call_worker(args):
    region_obs, pdict, tau, n_min = args
    return call_sites(region_obs, HMMParams.from_dict(pdict), tau=tau, n_min=n_min)


def parallel_call(
    region_obs: Sequence[Tuple[Region, TranscriptObs]],
    params: HMMParams,
    tau: float = 0.8,
    n_min: int = 10,
    threads: int = 1,
) -> pd.DataFrame:
    if threads <= 1 or len(region_obs) < 2:
        return call_sites(region_obs, params, tau=tau, n_min=n_min)
    jobs = [(c, params.to_dict(), tau, n_min) for c in _chunks(list(region_obs), threads * 4)]
    with ProcessPoolExecutor(max_workers=threads) as ex:
        parts = [df for df in ex.map(_call_worker, jobs) if len(df)]
    if not parts:
        return call_sites([], params, tau=tau, n_min=n_min)
    return pd.concat(parts, ignore_index=True)


def training_subset(
    region_obs: Sequence[Tuple[Region, TranscriptObs]],
    max_nt: int,
    seed: int = 0,
) -> List[TranscriptObs]:
    """Random regions totalling about `max_nt` nucleotides (0 = all).

    Baum-Welch cost grows with total chain length, and the parameters are
    global, so a few Mb of covered sequence is plenty to fit them.
    """
    trs = [tr for _, tr in region_obs if tr.A_idx.size]
    if max_nt <= 0 or sum(tr.L for tr in trs) <= max_nt:
        return trs
    rng = np.random.default_rng(seed)
    picked, total = [], 0
    for k in rng.permutation(len(trs)):
        picked.append(trs[k])
        total += trs[k].L
        if total >= max_nt:
            break
    return picked


def default_threads() -> int:
    return max(1, min(8, os.cpu_count() or 1))
