"""Generate a small, reproducible GLORI-seq example in the GLORI-tools layout.

Writes, under this directory:

    full.clean.fq                      clean single-end reads (sense strand)
    truth.bed                          true m6A sites (validation only)
    ref/index/genome/syn.main.fa(.fai)             unconverted genome
    ref/index/genome/syn.AG_conversion.fa(.fai)    A->G genome   (+ STAR dir)
    ref/index/genome/syn.rvsCom.fa                 T->C genome   (+ STAR dir)
    ref/index/transcriptome/syn_rna2.fa.AG_conversion.fa  (+ bowtie *.ebwt)
    ref/anno/syn.gtf_change2Ens.tbl2
    ref/anno/syn.gtf_change2Ens.tbl2.noredundance.base

Five genes on three contigs, both strands, two of them spliced, so every
alignment stage (STAR on both converted genomes, bowtie on the
transcriptome) and the transcript-to-genome lift get exercised.

Per adenosine, reads keep A with a state-specific probability:
    unmethylated ~0.02, true m6A ~0.70 (mostly DRACH), artifact patch ~0.25.

    python build_example.py              # data + indices (needs STAR, bowtie)
    python build_example.py --no-index   # data only
"""
from __future__ import annotations

import argparse
import math
import random
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pysam

HERE = Path(__file__).resolve().parent
SEED = 20251121

CONTIGS = {"chr1": 3200, "chr2": 3100, "chr3": 3700}
# (transcript, gene, contig, strand, exons as 0-based half-open genomic intervals)
GENES = [
    ("NM_000001", "GENE1", "chr1", "+", [(300, 1300)]),
    ("NM_000002", "GENE2", "chr1", "-", [(2000, 2900)]),
    ("NM_000003", "GENE3", "chr2", "+", [(200, 500), (1300, 1600), (2400, 2800)]),
    ("NM_000004", "GENE4", "chr3", "-", [(400, 800), (1700, 2100)]),
    ("NM_000005", "GENE5", "chr3", "+", [(2600, 3400)]),
]
PARAMS = {"U": (0.02, 60.0), "M": (0.70, 30.0), "B": (0.25, 25.0)}

_RC = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def is_drach(seq: str, i: int) -> bool:
    if i < 2 or i > len(seq) - 3:
        return False
    a, b, c, d, e = seq[i - 2:i + 3]
    return a in "AGT" and b in "AG" and c == "A" and d == "C" and e in "ACT"


def write_fasta(path: Path, recs: dict) -> None:
    with open(path, "w") as f:
        for name, seq in recs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")


def make_genome(rng: random.Random):
    genome = {c: [rng.choice("ACGT") for _ in range(L)] for c, L in CONTIGS.items()}
    tx_seqs = {}
    for tx, gene, contig, strand, exons in GENES:
        L = sum(e - s for s, e in exons)
        seq = [rng.choices("ACGT", weights=[27, 23, 23, 27])[0] for _ in range(L)]
        for p in rng.sample(range(2, L - 3), max(3, L // 120)):
            seq[p - 2:p + 3] = list("GGACT")          # canonical DRACH
        tseq = "".join(seq)
        tx_seqs[tx] = tseq
        plus = tseq if strand == "+" else revcomp(tseq)
        off = 0
        for s, e in exons:
            genome[contig][s:e] = list(plus[off:off + e - s])
            off += e - s
        for (s0, e0), (s1, _) in zip(exons, exons[1:]):  # intron GT..AG on gene strand
            donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
            genome[contig][e0:e0 + 2] = list(donor)
            genome[contig][s1 - 2:s1] = list(acceptor)
    return {c: "".join(v) for c, v in genome.items()}, tx_seqs


def tx_to_genome(exons, strand, t):
    L = sum(e - s for s, e in exons)
    a = t if strand == "+" else L - 1 - t
    for s, e in exons:
        if a < e - s:
            return s + a
        a -= e - s
    raise IndexError(t)


def simulate(rng: random.Random, nprng: np.random.Generator, tx_seqs, coverage=40):
    reads, truth = [], []
    for tx, gene, contig, strand, exons in GENES:
        seq = tx_seqs[tx]
        L = len(seq)
        drach = [i for i in range(L) if seq[i] == "A" and is_drach(seq, i)]
        other = [i for i in range(L) if seq[i] == "A" and not is_drach(seq, i)]
        sites = set(rng.sample(drach, max(2, len(drach) // 3))) | set(rng.sample(other, 1))
        a0 = rng.randint(L // 4, L // 2)
        artifact = set(range(a0, min(L, a0 + rng.randint(60, 90))))
        sites -= artifact
        sites = {i for i in sites if 20 <= i < L - 20}
        pA = np.zeros(L)
        for i in range(L):
            if seq[i] == "A":
                st = "M" if i in sites else "B" if i in artifact else "U"
                mu, k = PARAMS[st]
                pA[i] = nprng.beta(mu * k, (1 - mu) * k)
        n_reads = int(coverage * L / 75)
        for _ in range(n_reads):
            rl = rng.randint(50, 100)
            st = rng.randint(0, L - rl)
            r = [("G" if (b == "A" and nprng.random() > pA[st + k]) else b)
                 for k, b in enumerate(seq[st:st + rl])]
            reads.append("".join(r))
        for i in sites:
            truth.append((contig, tx_to_genome(exons, strand, i), strand, gene))
    for _ in range(50):                                   # unalignable junk
        reads.append("".join(rng.choice("ACGT") for _ in range(rng.randint(50, 100))))
    rng.shuffle(reads)
    return reads, sorted(truth)


def write_dataset(out: Path = HERE) -> Path:
    rng, nprng = random.Random(SEED), np.random.default_rng(SEED)
    genome, tx_seqs = make_genome(rng)
    reads, truth = simulate(rng, nprng, tx_seqs)

    gdir, tdir, adir = out / "ref/index/genome", out / "ref/index/transcriptome", out / "ref/anno"
    for d in (gdir, tdir, adir):
        d.mkdir(parents=True, exist_ok=True)

    write_fasta(gdir / "syn.main.fa", genome)
    write_fasta(gdir / "syn.AG_conversion.fa",
                {f"{c}_AG_converted": s.replace("A", "G") for c, s in genome.items()})
    write_fasta(gdir / "syn.rvsCom.fa",
                {f"{c}_AG_converted": s.replace("T", "C") for c, s in genome.items()})
    pysam.faidx(str(gdir / "syn.main.fa"))
    pysam.faidx(str(gdir / "syn.AG_conversion.fa"))
    write_fasta(tdir / "syn_rna2.fa.AG_conversion.fa",
                {f"{tx}_AG_converted": tx_seqs[tx].replace("A", "G") for tx, *_ in GENES})

    with open(adir / "syn.gtf_change2Ens.tbl2", "w") as f:
        for tx, gene, contig, strand, exons in GENES:
            starts = "".join(f"{s}," for s, _ in exons)
            ends = "".join(f"{e}," for _, e in exons)
            f.write("\t".join(map(str, [
                0, tx, contig, strand, exons[0][0], exons[-1][1], exons[0][0], exons[-1][1],
                len(exons), starts, ends, 0, gene, "cmpl", "cmpl", "0,"])) + "\n")
    with open(adir / "syn.gtf_change2Ens.tbl2.noredundance.base", "w") as f:
        for tx, gene, contig, strand, exons in GENES:
            for s, e in exons:
                for p in range(s, e):
                    f.write(f"{contig}\t{p}\t{p + 1}\t{strand}\t{gene}\t{tx}\t{gene}\t{tx}\tprotein_coding\n")

    with open(out / "full.clean.fq", "w") as f:
        for k, r in enumerate(reads):
            f.write(f"@SIM.{k} length={len(r)}\n{r}\n+\n{'I' * len(r)}\n")
    with open(out / "truth.bed", "w") as f:
        for contig, pos, strand, gene in truth:
            f.write(f"{contig}\t{pos}\t{pos + 1}\t{gene}\t1000\t{strand}\n")
    print(f"wrote {len(reads)} reads, {len(truth)} true sites under {out}")
    return out / "ref"


def build_indices(ref: Path, threads: int = 2) -> None:
    gdir, tdir = ref / "index/genome", ref / "index/transcriptome"
    gsize = sum(CONTIGS.values())
    nbases = int(math.floor(min(14, math.log2(gsize) / 2 - 1)))  # STAR manual formula
    for name in ("syn.AG_conversion", "syn.rvsCom"):
        d = gdir / name
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir()
        subprocess.run(
            ["STAR", "--runMode", "genomeGenerate", "--runThreadN", str(threads),
             "--genomeDir", str(d), "--genomeFastaFiles", str(gdir / f"{name}.fa"),
             "--genomeSAindexNbases", str(nbases), "--outFileNamePrefix", str(d) + "/"],
            check=True, stdout=subprocess.DEVNULL)
    fa = tdir / "syn_rna2.fa.AG_conversion.fa"
    subprocess.run(["bowtie-build", "-q", str(fa), str(fa)], check=True)
    print(f"built STAR and bowtie indices under {ref / 'index'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-index", action="store_true")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    ref = write_dataset()
    if not a.no_index:
        build_indices(ref, a.threads)
