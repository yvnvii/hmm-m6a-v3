"""Command-line interface: clean GLORI-seq FASTQ + GLORI-tools references -> m6A calls.

Inputs (GLORI-tools layout; pass --ref-dir to find them automatically)
    full.clean.fq                                      cleaned reads
    index/genome/hg38.AG_conversion.fa (+ .fai, STAR dir hg38.AG_conversion/)
    index/genome/hg38.rvsCom.fa                       (+ STAR dir hg38.rvsCom/)
    index/genome/hg38.main.fa (+ .fai)                 unconverted genome
    index/transcriptome/<...>_rna2.fa.AG_conversion.fa (+ bowtie *.ebwt)
    anno/<...>.tbl2                                    transcript exon table
    anno/<...>.tbl2.noredundance.base                  per-base gene labels (optional)

Subcommands
    hmm-m6a-fastq run    --fastq full.clean.fq --ref-dir REF --outdir out/
    hmm-m6a-fastq run-bam --bam GLORI/x_merged.sorted.bam --ref-dir REF --outdir out/
    hmm-m6a-fastq align  --fastq full.clean.fq --ref-dir REF -o sample.bam
    hmm-m6a-fastq pileup --bam sample.bam --ref-dir REF -o sample.pileup.tsv
    hmm-m6a-fastq train  --pileup sample.pileup.tsv -o params.json
    hmm-m6a-fastq call   --pileup sample.pileup.tsv --params params.json --ref-dir REF -o calls.tsv
    hmm-m6a-fastq default-params --states {3,4} -o params.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from hmm_m6a_fastq import __version__
from hmm_m6a_fastq.align import AlignOptions, GloriRefs, align_fastq
from hmm_m6a_fastq.annotation import annotate_sites
from hmm_m6a_fastq.hmm import HMMParams
from hmm_m6a_fastq.pileup import Region, _has_bai, regions_from_bed
from hmm_m6a_fastq.regions import (
    GLORI_CONTIG_SUFFIX,
    default_threads,
    parallel_call,
    parallel_pileup,
    regions_from_coverage,
    regions_to_bed,
    training_subset,
    with_bam_contigs,
)
from hmm_m6a_fastq.train import TranscriptObs, baum_welch


# ---------------------------------------------------------------------------
# Pileup TSV (de)serialisation
# ---------------------------------------------------------------------------
def _pileup_to_tsv(
    region_obs: Sequence[Tuple[Region, TranscriptObs]],
    out_path: str,
) -> None:
    """One row per adenosine with all the data the HMM needs.

    Columns: region, contig, strand, pos0_in_region, n, A, drach, region_L,
             region_start, region_end
    """
    rows: List[dict] = []
    for region, tr in region_obs:
        for i in tr.A_idx:
            i = int(i)
            rows.append({
                "region": region.name,
                "contig": region.contig,
                "strand": region.strand,
                "region_start": region.start,
                "region_end": region.end,
                "region_L": tr.L,
                "pos0_in_region": i,
                "A": int(tr.A_count[i]),
                "n": int(tr.n[i]),
                "drach": int(bool(tr.drach[i])),
            })
    cols = ["region", "contig", "strand", "region_start", "region_end", "region_L",
            "pos0_in_region", "A", "n", "drach"]
    pd.DataFrame(rows, columns=cols).to_csv(out_path, sep="\t", index=False)


def _tsv_to_pileup(path: str) -> List[Tuple[Region, TranscriptObs]]:
    """Reverse of `_pileup_to_tsv`."""
    df = pd.read_csv(path, sep="\t", dtype={"contig": str})
    out: List[Tuple[Region, TranscriptObs]] = []
    for region_name, sub in df.groupby("region", sort=False):
        sub = sub.sort_values("pos0_in_region")
        contig = sub["contig"].iloc[0]
        strand = sub["strand"].iloc[0]
        L = int(sub["region_L"].iloc[0])
        start = int(sub["region_start"].iloc[0])
        end = int(sub["region_end"].iloc[0])
        A_idx = sub["pos0_in_region"].astype(np.int32).to_numpy()
        A_count = np.zeros(L, dtype=np.int32)
        n = np.zeros(L, dtype=np.int32)
        drach = np.zeros(L, dtype=bool)
        A_count[A_idx] = sub["A"].astype(np.int32).to_numpy()
        n[A_idx] = sub["n"].astype(np.int32).to_numpy()
        drach[A_idx] = sub["drach"].astype(bool).to_numpy()
        region = Region(name=str(region_name), contig=contig, start=start, end=end, strand=strand)
        out.append((region, TranscriptObs(L=L, A_idx=A_idx, A_count=A_count, n=n, drach=drach)))
    return out


# ---------------------------------------------------------------------------
# Shared steps
# ---------------------------------------------------------------------------
def _log(args: argparse.Namespace, msg: str) -> None:
    if not args.quiet:
        print(msg, file=sys.stderr)


def _find_one(ref_dir: Optional[str], pattern: str, flag: str, required: bool) -> Optional[str]:
    if ref_dir:
        hits = sorted(glob.glob(os.path.join(ref_dir, pattern)))
        if len(hits) == 1:
            return hits[0]
        if hits:
            raise SystemExit(f"several {pattern} under {ref_dir}; pass {flag}")
    if required:
        raise SystemExit(f"{flag} is required (or --ref-dir containing {pattern})")
    return None


def _genome(args) -> str:
    return args.genome or _find_one(args.ref_dir, "index/genome/*.main.fa", "--genome", True)


def _baseanno(args) -> Optional[str]:
    return args.baseanno or _find_one(
        args.ref_dir, "anno/*.tbl2.noredundance.base", "--baseanno", False)


def _refs(args) -> GloriRefs:
    try:
        return GloriRefs.resolve(
            args.ref_dir, genome=args.genome, genome_ag=args.genome_ag,
            genome_rvs=args.genome_rvs, transcriptome=args.transcriptome,
            anno=args.anno, baseanno=args.baseanno)
    except (ValueError, FileNotFoundError) as e:
        raise SystemExit(str(e))


def _align_options(args) -> AlignOptions:
    return AlignOptions(
        threads=args.threads, mismatches=args.mismatches, min_frac=args.min_frac,
        star_extra=args.star_args.split(), bowtie_extra=args.bowtie_args.split(),
        glori_compat=args.glori_compat)


# Pileup filter defaults: the method's own, and GLORI-tools' (--glori-compat).
# A value given on the command line always wins over either.
PILEUP_DEFAULTS = {"min_mapq": 20, "min_baseq": 20, "max_unconverted": -1}
GLORI_PILEUP_DEFAULTS = {"min_mapq": 0, "min_baseq": 10, "max_unconverted": 3}


def _pileup_filters(args) -> dict:
    base = GLORI_PILEUP_DEFAULTS if args.glori_compat else PILEUP_DEFAULTS
    return {k: (getattr(args, k) if getattr(args, k) is not None else v)
            for k, v in base.items()}


def _do_pileup(args, bam: str, genome: str, regions_out: Optional[str] = None,
               strip_suffix: Optional[str] = None):
    """Regions -> per-A counts. `strip_suffix` is removed from the BAM's contig
    names to match the genome FASTA and in all outputs (the BAM is only read)."""
    filt = _pileup_filters(args)
    _log(args, "[pileup] filters: " + ", ".join(f"{k}={v}" for k, v in filt.items())
         + (" (GLORI-tools compatible)" if args.glori_compat else ""))
    if args.bed:
        regions = regions_from_bed(args.bed)
        if strip_suffix:
            regions = with_bam_contigs(regions, bam, strip_suffix)
    else:
        regions = regions_from_coverage(
            bam, min_reads=args.min_reads, max_gap=args.max_gap, min_mapq=filt["min_mapq"],
            strip_suffix=strip_suffix)
    total = sum(r.end - r.start for r in regions)
    _log(args, f"[pileup] {len(regions)} regions, {total:,} nt")
    if regions_out:
        regions_to_bed(regions, regions_out)
    return parallel_pileup(
        bam, genome, regions, threads=args.threads,
        glori_compat=args.glori_compat, **filt)


def _do_train(args, region_obs):
    transcripts = training_subset(region_obs, args.train_max_nt, seed=args.seed)
    _log(args, f"[train] Baum-Welch on {len(transcripts)} regions, "
               f"{sum(t.L for t in transcripts):,} nt")
    return baum_welch(
        transcripts, HMMParams.default(n_states=args.states),
        n_iter=args.n_iter, lam=args.lam, tol=args.tol, verbose=not args.quiet)


def _write_params(path, params: HMMParams, history=None) -> None:
    out = {"params": params.to_dict(), "version": __version__}
    if history is not None:
        out["history"] = {"logZ": list(map(float, history.logZ))}
    Path(path).write_text(json.dumps(out, indent=2))


def _load_params(path: str) -> HMMParams:
    return HMMParams.from_dict(json.loads(Path(path).read_text())["params"])


def _calls_to_bed(df: pd.DataFrame) -> str:
    called = df[df["call"]] if len(df) else df
    lines = []
    for r in called.itertuples(index=False):
        score = max(0, min(1000, int(round(1000.0 * float(r.pM)))))
        name = f"m6A;gene={r.gene};n={int(r.n)};A={int(r.A)};drach={int(bool(r.drach))}"
        lines.append(f"{r.contig}\t{int(r.pos0)}\t{int(r.pos0) + 1}\t{name}\t{score}\t{r.strand}")
    return "\n".join(lines) + ("\n" if lines else "")


def _do_call(args, region_obs, params: HMMParams, tsv, bed, baseanno: Optional[str]) -> None:
    df = parallel_call(region_obs, params, tau=args.tau, n_min=args.n_min, threads=args.threads)
    df = annotate_sites(df, baseanno, quiet=args.quiet)
    if len(df):
        df = df.sort_values(["contig", "pos0", "strand"], kind="stable")
    df.to_csv(tsv, sep="\t", index=False)
    n_called = int(df["call"].sum()) if len(df) else 0
    print(f"wrote {tsv} ({n_called} sites called)", file=sys.stderr)
    if bed:
        Path(bed).write_text(_calls_to_bed(df))
        print(f"wrote {bed}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------
def cmd_align(args) -> int:
    align_fastq(args.fastq, _refs(args), args.output, _align_options(args),
                tmp_dir=args.tmp_dir, keep_tmp=args.keep_tmp, quiet=args.quiet)
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


def cmd_pileup(args) -> int:
    region_obs = _do_pileup(args, args.bam, _genome(args), args.regions_out)
    _pileup_to_tsv(region_obs, args.output)
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


def cmd_train(args) -> int:
    params, history = _do_train(args, _tsv_to_pileup(args.pileup))
    _write_params(args.output, params, history)
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


def cmd_call(args) -> int:
    _do_call(args, _tsv_to_pileup(args.pileup), _load_params(args.params),
             args.output, args.bed_out, _baseanno(args))
    return 0


def cmd_default_params(args) -> int:
    _write_params(args.output, HMMParams.default(n_states=args.states))
    print(f"wrote {args.output} ({args.states}-state default parameters)", file=sys.stderr)
    return 0


def cmd_run(args) -> int:
    """Clean FASTQ -> BAM -> pileup -> [train] -> annotated calls."""
    refs = _refs(args)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    prefix = str(out / args.prefix)

    bam = f"{prefix}.bam"
    align_fastq(args.fastq, refs, bam, _align_options(args),
                tmp_dir=args.tmp_dir, keep_tmp=args.keep_tmp, quiet=args.quiet)
    return _score_bam(args, bam, refs.genome, refs.baseanno, prefix)


def _bam_prefix(bam: str) -> str:
    name = os.path.basename(bam)
    for suffix in ("_merged.sorted.bam", ".sorted.bam", ".bam"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def cmd_run_bam(args) -> int:
    """GLORI-tools BAM -> pileup -> [train] -> annotated calls (no alignment)."""
    if not os.path.exists(args.bam):
        raise SystemExit(f"no such BAM: {args.bam}")
    if not _has_bai(args.bam):
        import pysam
        _log(args, f"[run-bam] indexing {args.bam}")
        pysam.index(args.bam)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    prefix = str(out / (args.prefix or _bam_prefix(args.bam)))
    return _score_bam(args, args.bam, _genome(args), _baseanno(args), prefix,
                      strip_suffix=args.contig_suffix or None)


def _score_bam(args, bam: str, genome: str, baseanno: Optional[str], prefix: str,
               strip_suffix: Optional[str] = None) -> int:
    """Shared tail of `run` and `run-bam`: regions, pileup, train/params, calls."""
    region_obs = _do_pileup(args, bam, genome, f"{prefix}.regions.bed",
                            strip_suffix=strip_suffix)
    _pileup_to_tsv(region_obs, f"{prefix}.pileup.tsv")
    if not region_obs:
        print("no region reached --min-reads; nothing to score.", file=sys.stderr)
        return 1

    if args.params:
        params = _load_params(args.params)
        _log(args, f"loaded pretrained parameters from {args.params}")
    else:
        params, history = _do_train(args, region_obs)
        _write_params(f"{prefix}.params.json", params, history)

    _do_call(args, region_obs, params, f"{prefix}.calls.tsv", f"{prefix}.calls.bed",
             baseanno)
    return 0


# ---------------------------------------------------------------------------
# Argparse
# ---------------------------------------------------------------------------
def _add_ref_dir(p) -> None:
    p.add_argument("--ref-dir", default=None,
                   help="GLORI-tools reference root holding index/ and anno/; "
                        "any file not given explicitly is found there")


def _add_ref_args(p) -> None:
    _add_ref_dir(p)
    g = p.add_argument_group("reference files (override --ref-dir)")
    g.add_argument("--genome", help="index/genome/hg38.main.fa (unconverted, with .fai)")
    g.add_argument("--genome-ag", help="index/genome/hg38.AG_conversion.fa "
                                       "(STAR index in the folder of the same name)")
    g.add_argument("--genome-rvs", help="index/genome/hg38.rvsCom.fa "
                                        "(STAR index in the folder of the same name)")
    g.add_argument("--transcriptome", help="index/transcriptome/..._rna2.fa.AG_conversion.fa "
                                           "(bowtie index with the same prefix)")
    g.add_argument("--anno", help="anno/..._change2Ens.tbl2")
    g.add_argument("--baseanno", help="anno/..._change2Ens.tbl2.noredundance.base (optional)")


def _add_align_args(p) -> None:
    p.add_argument("--fastq", required=True, nargs="+",
                   help="clean single-end FASTQ (adapter/quality trimmed, deduplicated, "
                        "UMI removed); several files are concatenated")
    g = p.add_argument_group("alignment")
    g.add_argument("--mismatches", type=int, default=2,
                   help="max mismatches (STAR --outFilterMismatchNmax, bowtie -v)")
    g.add_argument("--min-frac", type=float, default=0.5,
                   help="STAR --outFilterScoreMinOverLread / --outFilterMatchNminOverLread")
    g.add_argument("--star-args", default="", help='extra STAR options, quoted')
    g.add_argument("--bowtie-args", default="", help='extra bowtie options, quoted')
    g.add_argument("--tmp-dir", default=None, help="scratch space for converted reads")
    g.add_argument("--keep-tmp", action="store_true", help="keep per-stage intermediates")


def _add_pileup_args(p) -> None:
    g = p.add_argument_group("pileup")
    g.add_argument("--bed", default=None,
                   help="optional BED6 of regions to score (one HMM chain per row); "
                        "default: build regions from read coverage")
    g.add_argument("--min-reads", type=int, default=10,
                   help="drop coverage regions with fewer reads (default 10)")
    g.add_argument("--max-gap", type=int, default=50,
                   help="merge reads closer than this into one region (default 50 nt)")
    g.add_argument("--min-mapq", type=int, default=None,
                   help="minimum MAPQ (default 20; --glori-compat: 0)")
    g.add_argument("--min-baseq", type=int, default=None,
                   help="minimum base quality at the site (default 20; --glori-compat: 10)")
    g.add_argument("--max-unconverted", type=int, default=None,
                   help="drop reads with more than this many unconverted A's; -1 = off "
                        "(default -1; --glori-compat: 3, counting every A in the read)")


def _add_compat_arg(p) -> None:
    p.add_argument("--glori-compat", action="store_true",
                   help="reproduce GLORI-tools' reads and per-site counts exactly, for a "
                        "like-for-like comparison with its binomial test: its alignment "
                        "strand rules and base restoration (including its bug on rvsCom "
                        "forward hits), whole-read A counting, the >= 2 reads rule, and "
                        "its filter values. Use it for both `align` and `pileup`.")


def _add_train_args(p) -> None:
    g = p.add_argument_group("training")
    g.add_argument("--states", type=int, default=3, choices=[3, 4],
                   help="3 = U/M/B (default); 4 adds N, a non-convertible state")
    g.add_argument("--n-iter", type=int, default=30)
    g.add_argument("--lam", type=float, default=1.0, help="Dirichlet pseudocount")
    g.add_argument("--tol", type=float, default=1e-3,
                   help="relative log-likelihood convergence tolerance")
    g.add_argument("--train-max-nt", type=int, default=2_000_000,
                   help="train on random regions totalling this many nt (0 = all)")
    g.add_argument("--seed", type=int, default=0, help="seed for the training subset")


def _add_call_args(p) -> None:
    g = p.add_argument_group("calling")
    g.add_argument("--tau", type=float, default=0.8, help="posterior threshold for M")
    g.add_argument("--n-min", type=int, default=10, help="minimum informative coverage")


def _add_threads(p) -> None:
    p.add_argument("--threads", type=int, default=default_threads())


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hmm-m6a-fastq",
        description="Per-nucleotide HMM for m6A detection from clean GLORI-seq FASTQ.",
    )
    p.add_argument("--version", action="version", version=f"hmm-m6a-fastq {__version__}")
    p.add_argument("-q", "--quiet", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    sr = sub.add_parser("run", help="end-to-end: clean FASTQ -> annotated m6A calls")
    _add_align_args(sr)
    _add_ref_args(sr)
    _add_pileup_args(sr)
    _add_compat_arg(sr)
    _add_train_args(sr)
    _add_call_args(sr)
    _add_threads(sr)
    sr.add_argument("--params", default=None,
                    help="JSON of pretrained parameters; skips training if given")
    sr.add_argument("--outdir", required=True)
    sr.add_argument("--prefix", default="sample", help="output file prefix")
    sr.set_defaults(func=cmd_run)

    sb = sub.add_parser(
        "run-bam", help="end-to-end from a GLORI-tools BAM (no alignment) -> annotated calls",
        description="Start from GLORI-tools' <prefix>_merged.sorted.bam: regions -> pileup "
                    "-> train (or --params) -> call -> gene labels, with the same outputs as "
                    "`run`. GLORI-tools' counting rules (--glori-compat) are on by default; "
                    "explicit filter options still override them.")
    sb.add_argument("--bam", required=True,
                    help="GLORI-tools <prefix>_merged.sorted.bam (indexed if no .bai)")
    _add_ref_dir(sb)
    sb.add_argument("--genome", help="index/genome/hg38.main.fa (unconverted, with .fai)")
    sb.add_argument("--baseanno", help="anno/..._change2Ens.tbl2.noredundance.base (optional)")
    sb.add_argument("--contig-suffix", default=GLORI_CONTIG_SUFFIX,
                    help=f"suffix stripped from BAM contig names to match the genome and in "
                         f"all outputs (default {GLORI_CONTIG_SUFFIX}; '' = none)")
    _add_pileup_args(sb)
    sb.add_argument("--glori-compat", action=argparse.BooleanOptionalAction, default=True,
                    help="GLORI-tools' counting rules: whole-read A count, >= 2 reads per "
                         "position, mapq 0 / baseq 10 / max-unconverted 3 (default on)")
    _add_train_args(sb)
    _add_call_args(sb)
    _add_threads(sb)
    sb.add_argument("--params", default=None,
                    help="JSON of pretrained parameters; skips training if given")
    sb.add_argument("--outdir", required=True)
    sb.add_argument("--prefix", default=None,
                    help="output file prefix (default: the BAM name without "
                         "_merged.sorted.bam)")
    sb.set_defaults(func=cmd_run_bam)

    sa = sub.add_parser("align", help="clean FASTQ -> genome BAM (STAR, STAR, bowtie)")
    _add_align_args(sa)
    _add_ref_args(sa)
    _add_compat_arg(sa)
    _add_threads(sa)
    sa.add_argument("-o", "--output", required=True, help="output BAM path")
    sa.set_defaults(func=cmd_align)

    sp = sub.add_parser("pileup", help="BAM (from `align`) -> per-A pileup TSV")
    sp.add_argument("--bam", required=True)
    _add_ref_dir(sp)
    sp.add_argument("--genome", help="index/genome/hg38.main.fa")
    _add_pileup_args(sp)
    _add_compat_arg(sp)
    _add_threads(sp)
    sp.add_argument("--regions-out", default=None, help="also write the regions as BED6")
    sp.add_argument("-o", "--output", required=True, help="output TSV path")
    sp.set_defaults(func=cmd_pileup)

    st = sub.add_parser("train", help="pileup TSV -> Baum-Welch JSON parameters")
    st.add_argument("--pileup", required=True)
    _add_train_args(st)
    st.add_argument("-o", "--output", required=True, help="output JSON path")
    st.set_defaults(func=cmd_train)

    sc = sub.add_parser("call", help="pileup TSV + params -> annotated per-site calls")
    sc.add_argument("--pileup", required=True)
    sc.add_argument("--params", required=True)
    _add_ref_dir(sc)
    sc.add_argument("--baseanno", help="anno/...tbl2.noredundance.base (optional)")
    _add_call_args(sc)
    _add_threads(sc)
    sc.add_argument("-o", "--output", required=True, help="output TSV path")
    sc.add_argument("--bed-out", default=None, help="optional BED6 of called sites")
    sc.set_defaults(func=cmd_call)

    sd = sub.add_parser("default-params", help="write default (untrained) parameters to JSON")
    sd.add_argument("--states", type=int, default=3, choices=[3, 4])
    sd.add_argument("-o", "--output", required=True, help="output JSON path")
    sd.set_defaults(func=cmd_default_params)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
