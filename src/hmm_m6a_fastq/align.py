"""Conversion-aware alignment of clean GLORI-seq reads against the GLORI-tools
reference layout (same cascade as GLORI-tools' mapping_reads.py).

Reads are sense-strand single-end GLORI reads (A-less library). Every A in a
read is changed to G in silico, then:

    stage 1  STAR   -> hg38.AG_conversion   (A->G genome)
             keep forward-strand hits        = reads from + strand genes
    stage 2  STAR   -> hg38.rvsCom          (T->C genome, forward coords)
             reads not kept in stage 1; keep reverse-strand hits
                                             = reads from - strand genes
    stage 3  bowtie -> transcriptome AG_conversion (longest transcript/gene)
             reads not kept in stages 1-2; keep forward hits, lift them to
             the genome with the tbl2 exon table (spliced reads end up here)

Unique hits only (STAR --outFilterMultimapNmax 1, bowtie -m 1). A read that
is not kept at one stage (unmapped, too many loci, or wrong strand) is tried
at the next. With `AlignOptions.glori_compat` the strand rules and base
restoration are GLORI-tools' instead (see `_process_stage`), so the BAM holds
exactly the reads GLORI-tools keeps.

Kept reads get their ORIGINAL bases back so the pileup sees real A
(retained) vs G (converted) calls, and contig names lose the
"_AG_converted" suffix so they match hg38.main.fa.

Restoring bases needs no lookup table: each converted read is renamed
"<ordinal>_<hex mask of its original A positions>", which travels through
STAR and bowtie untouched.
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import IO, Dict, Iterator, List, Optional, Sequence, Tuple

import pysam

from hmm_m6a_fastq.annotation import TX_SUFFIX, TranscriptModel, load_tbl2

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_TO_G = str.maketrans("A", "G")


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise RuntimeError(
            f"'{tool}' not found on PATH (needs STAR >= 2.7.5c and bowtie 1.x, "
            f"e.g. `conda install -c bioconda star bowtie`).")
    return path


# ---------------------------------------------------------------------------
# Reference layout
# ---------------------------------------------------------------------------
@dataclass
class GloriRefs:
    """Paths of the GLORI-tools reference files.

    genome         index/genome/<g>.main.fa           unconverted (+ .fai)
    genome_ag      index/genome/<g>.AG_conversion.fa  STAR index = same path minus ".fa"
    genome_rvs     index/genome/<g>.rvsCom.fa         STAR index = same path minus ".fa"
    transcriptome  index/transcriptome/<t>.AG_conversion.fa   bowtie index prefix
    anno           anno/<a>.tbl2                       transcript exon table
    baseanno       anno/<a>.tbl2.noredundance.base     per-base gene labels (optional)
    """

    genome: str
    genome_ag: str
    genome_rvs: str
    transcriptome: str
    anno: str
    baseanno: Optional[str] = None

    @staticmethod
    def resolve(ref_dir: Optional[str] = None, **explicit: Optional[str]) -> "GloriRefs":
        """Fill unspecified paths from the standard layout under `ref_dir`."""

        def find(pattern: str, what: str, required: bool = True) -> Optional[str]:
            if ref_dir is None:
                if required:
                    raise ValueError(f"--{what.replace('_', '-')} is required (or pass --ref-dir)")
                return None
            hits = sorted(glob.glob(os.path.join(ref_dir, pattern)))
            if len(hits) == 1:
                return hits[0]
            if not hits and not required:
                return None
            raise ValueError(
                f"expected exactly one {pattern} under {ref_dir} for {what}, "
                f"found {len(hits)}; pass --{what.replace('_', '-')} explicitly")

        def pick(key: str, pattern: str, required: bool = True) -> Optional[str]:
            return explicit.get(key) or find(pattern, key, required)

        refs = GloriRefs(
            genome=pick("genome", "index/genome/*.main.fa"),
            genome_ag=pick("genome_ag", "index/genome/*.AG_conversion.fa"),
            genome_rvs=pick("genome_rvs", "index/genome/*.rvsCom.fa"),
            transcriptome=pick("transcriptome", "index/transcriptome/*.AG_conversion.fa"),
            anno=pick("anno", "anno/*.tbl2"),
            baseanno=pick("baseanno", "anno/*.tbl2.noredundance.base", required=False),
        )
        refs.check()
        return refs

    @staticmethod
    def star_dir(fa: str) -> str:
        return fa[:-3] if fa.endswith(".fa") else fa

    def check(self) -> None:
        problems = []
        for label, p in (("genome", self.genome), ("anno", self.anno)):
            if not os.path.isfile(p):
                problems.append(f"{label}: {p} not found")
        if not os.path.isfile(self.genome + ".fai"):
            problems.append(f"genome index {self.genome}.fai not found (samtools faidx)")
        for label, fa in (("genome_ag", self.genome_ag), ("genome_rvs", self.genome_rvs)):
            d = self.star_dir(fa)
            if not os.path.isfile(os.path.join(d, "SA")):
                problems.append(f"{label}: no STAR index in {d}/")
        if not (glob.glob(self.transcriptome + ".1.ebwt") or glob.glob(self.transcriptome + ".1.ebwtl")):
            problems.append(f"transcriptome: no bowtie index {self.transcriptome}.1.ebwt")
        if self.baseanno and not os.path.isfile(self.baseanno):
            problems.append(f"baseanno: {self.baseanno} not found")
        if problems:
            raise FileNotFoundError("reference files missing:\n  " + "\n  ".join(problems))


# ---------------------------------------------------------------------------
# FASTQ I/O and read encoding
# ---------------------------------------------------------------------------
def _open_text(path: str) -> IO[str]:
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def iter_fastq(paths: Sequence[str]) -> Iterator[Tuple[str, str, str]]:
    """Yield (name, sequence, quality) from one or more FASTQ(.gz) files, in order."""
    for path in paths:
        with _open_text(path) as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline().rstrip("\r\n")
                fh.readline()
                q = fh.readline().rstrip("\r\n")
                if not h.startswith("@") or len(s) != len(q):
                    raise ValueError(f"malformed FASTQ record in {path}: {h.strip()!r}")
                yield h[1:].split()[0], s.upper(), q


def encode_name(i: int, seq: str) -> str:
    """Ordinal + hex bitmask of A positions (bit k = base k of the read)."""
    mask = 0
    for k, b in enumerate(seq):
        if b == "A":
            mask |= 1 << k
    return f"{i}_{mask:x}"


def restore_bases(conv_read: str, name: str) -> str:
    """Put the A's back into a converted read (read orientation)."""
    mask = int(name.rsplit("_", 1)[1], 16)
    if not mask:
        return conv_read
    s = list(conv_read)
    k = 0
    while mask:
        if mask & 1:
            s[k] = "A"
        mask >>= 1
        k += 1
    return "".join(s)


def restore_bases_mirrored(conv_read: str, name: str) -> str:
    """GLORI-tools' restoration of a stage-2 (rvsCom) forward-strand hit.

    mapping_reads.py reverseReads2 assumes every rvsCom hit is reverse: it
    reverse-complements the stored sequence, writes A at the read's original
    A positions and reverse-complements back. For a forward hit this puts T at
    the mirrored positions (len - 1 - k) instead of A at k. Reproduced as is so
    that both pipelines count the same bases.
    """
    mask = int(name.rsplit("_", 1)[1], 16)
    s = list(conv_read)
    L = len(s)
    k = 0
    while mask:
        if mask & 1:
            s[L - 1 - k] = "T"
        mask >>= 1
        k += 1
    return "".join(s)


def _write_fq(fh: IO[str], name: str, seq: str, qual: str) -> None:
    fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


# ---------------------------------------------------------------------------
# Aligners
# ---------------------------------------------------------------------------
@dataclass
class AlignOptions:
    threads: int = 4
    mismatches: int = 2          # STAR --outFilterMismatchNmax / bowtie -v
    min_frac: float = 0.5        # STAR --outFilterScoreMinOverLread / MatchNminOverLread
    star_extra: Sequence[str] = ()
    bowtie_extra: Sequence[str] = ()
    glori_compat: bool = False   # GLORI-tools strand rules and base restoration


def _run_star(fq: Path, star_dir: str, prefix: Path, opt: AlignOptions) -> Path:
    cmd = [
        _require("STAR"), "--runThreadN", str(opt.threads), "--genomeDir", star_dir,
        "--readFilesIn", str(fq), "--outFileNamePrefix", str(prefix),
        "--limitOutSJcollapsed", "5000000",
        "--outFilterMismatchNmax", str(opt.mismatches),
        "--outFilterScoreMinOverLread", str(opt.min_frac),
        "--outFilterMatchNminOverLread", str(opt.min_frac),
        "--seedSearchStartLmax", "30",
        "--outSAMattributes", "All", "--outSAMprimaryFlag", "AllBestScore",
        "--outMultimapperOrder", "Random", "--runRNGseed", "777",
        "--outSAMmultNmax", "1", "--outFilterMultimapNmax", "1",
        "--outSAMtype", "BAM", "Unsorted", "--outSAMunmapped", "Within",
        *opt.star_extra,
    ]
    with open(f"{prefix}stdout.log", "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    return Path(f"{prefix}Aligned.out.bam")


def _run_bowtie(fq: Path, index: str, sam: Path, opt: AlignOptions) -> Path:
    cmd = [
        _require("bowtie"), "-k", "1", "-m", "1", "-v", str(opt.mismatches),
        "--best", "--strata", "-p", str(opt.threads), "--sam",
        *opt.bowtie_extra, "-x", index, str(fq), str(sam),
    ]
    with open(f"{sam}.log", "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    return sam


# ---------------------------------------------------------------------------
# Record conversion
# ---------------------------------------------------------------------------
# MD/NM etc. describe the converted alignment, not the real one; keep only these.
_KEEP_TAGS = {"NH", "HI", "AS", "nM", "XM"}


def _read_orient(rec: pysam.AlignedSegment) -> Tuple[str, str]:
    """(sequence, quality string) in original read orientation."""
    seq = rec.query_sequence or ""
    qual = pysam.qualities_to_qualitystring(rec.query_qualities) if rec.query_qualities is not None \
        else "I" * len(seq)
    if rec.is_reverse:
        return revcomp(seq), qual[::-1]
    return seq, qual


def _emit(
    out_header: pysam.AlignmentHeader,
    tid: int,
    start: int,
    cigar: List[Tuple[int, int]],
    reverse: bool,
    read_seq: str,
    read_qual: str,
    name: str,
    mapq: int,
    tags: List[Tuple[str, object]],
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment(out_header)
    a.query_name = name
    a.flag = 16 if reverse else 0
    a.reference_id = tid
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigartuples = cigar
    if reverse:
        a.query_sequence = revcomp(read_seq)
        a.query_qualities = pysam.qualitystring_to_array(read_qual[::-1])
    else:
        a.query_sequence = read_seq
        a.query_qualities = pysam.qualitystring_to_array(read_qual)
    a.set_tags(tags)
    return a


@dataclass
class AlignStats:
    reads_in: int = 0
    genome_plus: int = 0          # stage 1, + strand genes
    genome_minus: int = 0         # stage 2, - strand genes (glori_compat: plus
                                  # GLORI-tools' rare forward hits)
    transcriptome: int = 0        # stage 3, lifted to genome
    wrong_strand: int = 0         # glori_compat only: stage 1/3 reverse hits dropped
    unaligned: int = 0            # no unique, accepted hit at any stage
    unlifted: int = 0             # stage-3 hits that could not be lifted
    off_genome: int = 0           # transcriptome hits on chroms missing from the main genome
    stage_logs: Dict[str, str] = field(default_factory=dict)

    @property
    def aligned(self) -> int:
        return self.genome_plus + self.genome_minus + self.transcriptome


def _process_stage(
    in_path: Path,
    stage: str,
    out_bam: pysam.AlignmentFile,
    next_fq: Optional[IO[str]],
    tids: Dict[str, int],
    stats: AlignStats,
    tx_models: Optional[Dict[str, TranscriptModel]] = None,
    glori_compat: bool = False,
) -> None:
    """Keep this stage's accepted reads (restored) and pass the rest on.

    Default: keep forward hits in ag and tx and reverse hits in rvs; every
    other read (unmapped or wrong strand) is retried at the next stage.

    glori_compat reproduces GLORI-tools' mapping_reads.py instead:
      ag  (STAR, A->G genome)  keep forward hits (samtools -F 20); reverse
                               hits are dropped, NOT retried, because only
                               STAR's unmapped reads go on to stage 2.
      rvs (STAR, T->C genome)  keep hits in either orientation (-F 4 -q 255);
                               forward hits get GLORI-tools' mirrored base
                               restoration (`restore_bases_mirrored`).
      tx  (bowtie)             keep forward hits (-F 20).
    """
    with pysam.AlignmentFile(str(in_path), "r", check_sq=False) as fh:
        for rec in fh:
            if rec.is_secondary or rec.is_supplementary:
                continue
            seq, qual = _read_orient(rec)
            name = rec.query_name
            kept = False
            mapped = not rec.is_unmapped
            if glori_compat:
                if mapped and rec.is_reverse and stage in ("ag", "tx"):
                    stats.wrong_strand += 1
                    continue
                accept = mapped
            else:
                accept = mapped and rec.is_reverse == (stage == "rvs")
            if accept:
                if glori_compat and stage == "rvs" and not rec.is_reverse:
                    read_seq = restore_bases_mirrored(seq, name)
                else:
                    read_seq = restore_bases(seq, name)
                ref = rec.reference_name
                base_ref = ref[: -len(TX_SUFFIX)] if ref.endswith(TX_SUFFIX) else ref
                tags = [(t, v) for t, v in rec.get_tags() if t in _KEEP_TAGS]
                if stage == "tx":
                    model = tx_models.get(base_ref) if tx_models else None
                    lifted = None
                    if model is not None and all(op == 0 for op, _ in rec.cigartuples):
                        lifted = model.lift(rec.reference_start, rec.query_length)
                    if lifted is None:
                        stats.unlifted += 1
                    elif model.chrom not in tids:
                        stats.off_genome += 1
                    else:
                        start, cigar = lifted
                        out_bam.write(_emit(
                            out_bam.header, tids[model.chrom], start, cigar,
                            model.strand == "-", read_seq, qual, name, rec.mapping_quality,
                            tags + [("TN", base_ref), ("GN", model.gene)]))
                        stats.transcriptome += 1
                        kept = True
                elif base_ref in tids:   # else: retried at the next stage
                    out_bam.write(_emit(
                        out_bam.header, tids[base_ref], rec.reference_start, rec.cigartuples,
                        rec.is_reverse, read_seq, qual, name, rec.mapping_quality, tags))
                    if stage == "ag":
                        stats.genome_plus += 1
                    else:
                        stats.genome_minus += 1
                    kept = True
            if not kept and next_fq is not None:
                _write_fq(next_fq, name, seq, qual)


def align_fastq(
    fastqs: Sequence[str],
    refs: GloriRefs,
    out_bam: str,
    opt: AlignOptions = AlignOptions(),
    tmp_dir: Optional[str] = None,
    keep_tmp: bool = False,
    quiet: bool = False,
) -> AlignStats:
    """Clean FASTQ -> sorted, indexed BAM of original-base reads on hg38.main.

    Also writes <out_bam>.stats.json and the aligners' logs next to the BAM.
    """
    out_bam = str(out_bam)
    Path(out_bam).parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="glori_align_", dir=tmp_dir or Path(out_bam).parent))
    stats = AlignStats()

    def say(msg: str) -> None:
        if not quiet:
            print(f"[align] {msg}", file=sys.stderr)

    try:
        # Final header: contigs of the unconverted genome.
        with pysam.FastaFile(refs.genome) as fa:
            sq = list(zip(fa.references, fa.lengths))
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": n, "LN": l} for n, l in sq],
            "PG": [{"ID": "hmm-m6a-fastq", "PN": "hmm-m6a-fastq",
                    "CL": "STAR(AG) -> STAR(rvsCom) -> bowtie(transcriptome)"}],
        })
        tids = {n: i for i, (n, _) in enumerate(sq)}

        conv_fq = work / "stage1.fq"
        with open(conv_fq, "w") as fo:
            for i, (_, s, q) in enumerate(iter_fastq(fastqs)):
                _write_fq(fo, encode_name(i, s), s.translate(_TO_G), q)
                stats.reads_in += 1
        say(f"{stats.reads_in} reads; A->G converted")

        unsorted = work / "merged.unsorted.bam"
        with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
            say("stage 1: STAR -> A->G genome (+ strand genes)")
            s1 = _run_star(conv_fq, refs.star_dir(refs.genome_ag), work / "ag.", opt)
            fq2 = work / "stage2.fq"
            with open(fq2, "w") as nf:
                _process_stage(s1, "ag", bam, nf, tids, stats, glori_compat=opt.glori_compat)
            s1.unlink()
            say(f"  kept {stats.genome_plus}")

            say("stage 2: STAR -> T->C genome (- strand genes)")
            s2 = _run_star(fq2, refs.star_dir(refs.genome_rvs), work / "rvs.", opt)
            fq3 = work / "stage3.fq"
            with open(fq3, "w") as nf:
                _process_stage(s2, "rvs", bam, nf, tids, stats, glori_compat=opt.glori_compat)
            s2.unlink()
            say(f"  kept {stats.genome_minus}")

            say("stage 3: bowtie -> transcriptome (lifted to genome via tbl2)")
            models = load_tbl2(refs.anno)
            s3 = _run_bowtie(fq3, refs.transcriptome, work / "tx.sam", opt)
            _process_stage(s3, "tx", bam, None, tids, stats, tx_models=models,
                           glori_compat=opt.glori_compat)
            s3.unlink()
            say(f"  kept {stats.transcriptome} (unliftable {stats.unlifted})")

        stats.unaligned = (stats.reads_in - stats.aligned - stats.wrong_strand
                           - stats.unlifted - stats.off_genome)
        for tag, src in (("STAR_AG", work / "ag.Log.final.out"),
                         ("STAR_rvsCom", work / "rvs.Log.final.out"),
                         ("bowtie_transcriptome", work / "tx.sam.log")):
            if src.exists():
                dst = f"{out_bam}.{tag}.log"
                shutil.copy(src, dst)
                stats.stage_logs[tag] = dst

        pysam.sort("-@", str(max(1, opt.threads)), "-T", str(work / "sort"),
                   "-o", out_bam, str(unsorted))
        pysam.index(out_bam)
    finally:
        if keep_tmp:
            say(f"intermediates kept in {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)

    Path(f"{out_bam}.stats.json").write_text(
        json.dumps({**asdict(stats), "aligned": stats.aligned,
                    "glori_compat": opt.glori_compat}, indent=2))
    n = max(stats.reads_in, 1)
    say(f"kept {stats.aligned}/{stats.reads_in} ({100 * stats.aligned / n:.1f}%): "
        f"+genome {stats.genome_plus}, -genome {stats.genome_minus}, "
        f"transcriptome {stats.transcriptome}; wrong strand {stats.wrong_strand}, "
        f"unaligned {stats.unaligned}")
    return stats
