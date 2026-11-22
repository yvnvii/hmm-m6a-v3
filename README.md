# hmm-m6a-fastq

hmm-m6a for **clean GLORI-seq FASTQ** input, using the **GLORI-tools
reference files** (converted genomes with STAR indices, a bowtie
transcriptome index, and the GLORI annotation tables). Reads are aligned
with the same three-stage cascade as GLORI-tools. Sites are then scored with
the hmm-m6a per-nucleotide HMM instead of GLORI's binomial test. With
`--glori-compat` the reads and per-site counts are exactly GLORI-tools', so
only the calling method differs (see [Comparing with GLORI-tools](#comparing-with-glori-tools)).
`hmm.py`, `mle.py`, `train.py` and the `train` subcommand are unchanged
from hmm-m6a.

```
full.clean.fq
   │  A→G in silico
   ├─ STAR   → hg38.AG_conversion   keep forward hits        (+ strand genes)
   ├─ STAR   → hg38.rvsCom          keep reverse hits        (− strand genes)
   └─ bowtie → transcriptome AG     keep forward hits, lift to genome via tbl2
   ▼
sample.bam (original bases, hg38.main coordinates)
   ▼  coverage regions → pileup against hg38.main.fa
per-A counts ──train (Baum-Welch)──▶ params ──call──▶ sites ──label genes (noredundance.base)
```

## Inputs

These are the only inputs:

| file | role |
|---|---|
| `full.clean.fq` | reads: adapter/quality trimmed (trim_galore), deduplicated on UMI + sequence (seqkit), UMI removed (fastx_trimmer). Plain-text or gzipped single-end FASTQ |
| `index/genome/hg38.AG_conversion.fa` (+`.fai`) | genome with every A→G; `index/genome/hg38.AG_conversion/` is its STAR index |
| `index/genome/hg38.rvsCom.fa` | genome with every T→C (forward coordinates); `index/genome/hg38.rvsCom/` is its STAR index |
| `index/genome/hg38.main.fa` (+`.fai`) | unconverted genome: true reference base at each position |
| `index/transcriptome/…_rna2.fa.AG_conversion.fa` (+`*.ebwt`) | longest transcript per gene, A→G, with bowtie index |
| `anno/…_change2Ens.tbl2` | transcript/exon table used to lift transcriptome hits to the genome |
| `anno/…_change2Ens.tbl2.noredundance.base` | per-base gene labels (optional; without it every site is `ELSE`) |

Pass the folder holding `index/` and `anno/` as `--ref-dir`, and each file
is found by its suffix. Any file can also be given explicitly
(`--genome`, `--genome-ag`, `--genome-rvs`, `--transcriptome`, `--anno`,
`--baseanno`), which overrides the lookup.

No regions BED is needed: HMM chains are built from read coverage (see
step 3). You can still pass `--bed` to restrict scoring.

## Installation

```bash
conda create -n hmm-m6a-fastq -c conda-forge -c bioconda python=3.11 star bowtie
conda activate hmm-m6a-fastq
pip install -e .
```

Needs STAR ≥ 2.7.5c and bowtie 1.x on `PATH`. Use the same STAR version
that built the indices.

## Quick start

```bash
python example/build_example.py      # synthetic reference + reads (needs STAR, bowtie)

hmm-m6a-fastq run \
  --fastq   example/full.clean.fq \
  --ref-dir example/ref \
  --outdir  example_out

python example/validate.py example_out/sample.calls.tsv
```

This takes about 5 seconds. Expected result: all 36 true sites called,
0 false positives, every site labelled with its gene.

On real data:

```bash
hmm-m6a-fastq run --fastq full.clean.fq --ref-dir /path/to/GLORI_ref \
                  --outdir out --prefix sample1 --threads 16
```

## Pipeline, step by step

### 1. `align`: clean FASTQ → genome BAM

| | |
|---|---|
| **input** | `--fastq`, reference files (`--ref-dir`) |
| **output** | `sample.bam` (+`.bai`), `sample.bam.stats.json`, `sample.bam.STAR_AG.log`, `.STAR_rvsCom.log`, `.bowtie_transcriptome.log` |

1. Every A in every read is changed to G. The positions of the original A's
   are stored in the read name, so no lookup table is needed later.
2. **STAR → AG genome.** Reads with a unique forward hit are kept (+ strand genes).
3. **STAR → rvsCom genome.** Reads not kept so far; those with a unique
   reverse hit are kept (− strand genes).
4. **bowtie → transcriptome.** Reads still not kept; unique forward hits
   are lifted to genome coordinates with the tbl2 exon table. Reads that
   cross an exon junction become spliced alignments (`N` in the CIGAR).
5. Kept reads get their original bases back, and contig names lose the
   `_AG_converted` suffix so they match `hg38.main.fa`.

Alignment settings match GLORI-tools: up to 2 mismatches,
`--outFilterScoreMinOverLread 0.5`, `--outFilterMatchNminOverLread 0.5`,
unique hits only (STAR `--outFilterMultimapNmax 1`, bowtie `-m 1`).

A read that hits the wrong strand at one stage is tried at the next.
GLORI-tools differs here; `--glori-compat` follows it instead (see below).

`stats.json` reports reads kept at each stage, plus unaligned and
unliftable reads (and, with `--glori-compat`, dropped wrong-strand reads).

### 2. `pileup`: BAM → per-adenosine counts

| | |
|---|---|
| **input** | `--bam`, `--genome` (hg38.main.fa) or `--ref-dir` |
| **output** | pileup TSV (one row per transcript-strand A), optional `--regions-out` BED |

- **Regions (HMM chains).** On each strand, aligned read blocks within
  `--max-gap` (50 nt) of each other are merged into one region. Regions
  with fewer than `--min-reads` (10) reads are dropped. Introns split
  regions, and uncovered sequence is skipped, so run time scales with
  covered sequence, not genome size.
- **Counting.** For every reference A on the region's strand (looked up
  in `hg38.main.fa`), the tool counts reads showing A and reads showing G,
  and flags A's in a DRACH motif.
- **Filters:** `--min-mapq` (20), `--min-baseq` (20) and
  `--max-unconverted N` (default `-1`, off). `--max-unconverted` drops reads
  with more than N unconverted A's, counted at aligned positions where the
  reference is A and the base passes `--min-baseq`. Columns are read to a
  depth of 10,000,000 (10,000 on `GL…` contigs), as in GLORI-tools; pysam's
  default of 8,000 would silently drop reads at deep loci.

### 3. `train`: Baum-Welch (unchanged from hmm-m6a)

| | |
|---|---|
| **input** | `--pileup`, `--states 3|4`, `--train-max-nt` |
| **output** | params JSON |

The parameters are shared across the whole genome, so training uses a
random subset of regions totalling `--train-max-nt` nucleotides (default
2 Mb; set 0 to use everything).

### 4. `call`: posteriors, calling, gene labels

| | |
|---|---|
| **input** | `--pileup`, `--params`, `--baseanno` or `--ref-dir` |
| **output** | calls TSV (every scored A) and optional `--bed-out` BED6 of called sites |

A site is called when all three hold:
- `P(M) ≥ --tau` (0.8)
- `P(M) > P(B)` (and `> P(N)` with 4 states)
- `n ≥ --n-min` (10)

Each row is then matched on (chrom, strand, position) against
`noredundance.base`, which adds the columns `gene`, `gene_name`,
`transcript` and `biotype`. Sites outside annotated exons get
`gene = ELSE`, as in GLORI-tools.

TSV columns: `region contig pos0 strand n A rate drach pU pM pB [pN] call gene gene_name transcript biotype`.
BED name field: `m6A;gene=…;n=…;A=…;drach=…`; score = `round(1000·pM)`.

### `run`: all of the above

| | |
|---|---|
| **input** | `--fastq`, `--ref-dir` (or explicit files), `--outdir` |
| **output** in `--outdir` | `sample.bam` (+ stats and logs), `sample.regions.bed`, `sample.pileup.tsv`, `sample.params.json`, `sample.calls.tsv`, `sample.calls.bed` |

Pass `--params` to skip training and score with existing parameters.
`--threads` sets the thread count for STAR, bowtie, pileup and calling.

## Comparing with GLORI-tools

To test the HMM against GLORI-tools' binomial test on identical input, pass
`--glori-compat` to `run`, or to both `align` and `pileup`. The reads in the
BAM and the per-site `n` and `A` are then exactly the `AGcov` and `Acov` that
GLORI-tools' caller tests at its default `--cutoff 3`:

| | default | `--glori-compat` (GLORI-tools source) |
|---|---|---|
| stage-1 / stage-3 reverse hits | tried at the next stage | dropped (`mapping_reads.py`, `samtools -F 20`) |
| stage-2 forward hits | tried at the next stage | kept, with GLORI-tools' base restoration, which writes T at mirrored positions instead of A (`reverseReads2`) |
| `--min-mapq` | 20 | 0 (no MAPQ filter) |
| `--min-baseq` | 20 | 10 (`pileup_genome_multiprocessing.py -q`) |
| `--max-unconverted` | -1 (off) | 3 (`m6A_caller.py --cutoff 3`) |
| how A's per read are counted | aligned reference A's that pass `--min-baseq` | every A in the read, soft clips and mismatches included |
| positions seen by < 2 reads | kept | n = 0 (`get_referbase.py`) |

The stage-2 restoration is a GLORI-tools bug, reproduced only so that both
pipelines count the same bases; keep it out of normal runs. A value given on
the command line overrides the compat default, so
`--glori-compat --max-unconverted -1` gives GLORI-tools' reads with the
read filter off (GLORI's unfiltered `AGcov`/`Acov` columns).

To confirm that both pipelines saw the same data, run GLORI-tools and v3 on
the same FASTQ and reference, then compare. The script exits 0 only if the
reads and every per-site count agree:

```bash
python scripts/compare_glori.py \
    --glori-format GLORI_out/x.totalformat.txt --v3-pileup out/x.pileup.tsv \
    --glori-bam GLORI_out/x_merged.sorted.bam --v3-bam out/x.bam
# v3 pileup made with --max-unconverted -1: compare against GLORI's unfiltered columns
python scripts/compare_glori.py --acut all --glori-format ... --v3-pileup ...
```

GLORI's site filters (coverage ≥ 15, A ≥ 5, ratio ≥ 0.1, `var_ratio`,
signal ratio, gene conversion rate) belong to its calling method and are
not applied. `--n-min 15` restricts the HMM to GLORI's coverage cutoff.

## Several samples: train once, call many

```bash
hmm-m6a-fastq run --fastq ref.clean.fq --ref-dir $REF --outdir out/ref --prefix ref
hmm-m6a-fastq run --fastq B.clean.fq   --ref-dir $REF --outdir out/B   --prefix B \
                  --params out/ref/ref.params.json
```

## Notes and limitations

- **Single-end, sense-strand reads only.** This matches GLORI's A-less
  library and the reference layout, which only has an A→G transcriptome.
- **Read names.** The BAM uses `<ordinal>_<hex>` read names. These encode
  each read's A positions, not the original FASTQ names.
- **Junction edges.** As in GLORI-tools, STAR sometimes extends a read a few
  nt into an intron or soft-clips it instead of splicing. On the synthetic
  data this affects 3 of ~190k aligned bases, none at a scored A.
- **Run time on hg38.** Alignment is dominated by STAR. The HMM scales with
  covered sequence: with `--threads`, pileup and calling run in parallel
  across regions, and training is capped by `--train-max-nt`.

## Layout

```
hmm-m6a-fastq/
├── src/hmm_m6a_fastq/
│   ├── align.py       STAR(AG) → STAR(rvsCom) → bowtie(transcriptome) cascade
│   ├── annotation.py  tbl2 lift-over, noredundance.base gene labels
│   ├── regions.py     coverage regions, parallel pileup/calling, training subset
│   ├── cli.py         run / align / pileup / train / call / default-params
│   ├── pileup.py      BAM + FASTA → per-A counts and DRACH   (from hmm-m6a)
│   ├── calling.py     posterior thresholding                 (from hmm-m6a)
│   ├── hmm.py         beta-binom PMF, FB, Viterbi, posterior  (unchanged)
│   ├── mle.py         weighted MLE for (mu, kappa)            (unchanged)
│   └── train.py       Baum-Welch with state-ordering constraint (unchanged)
├── example/           build_example.py (synthetic GLORI-tools layout), validate.py
├── scripts/           compare_glori.py (reads/counts vs GLORI-tools)
└── tests/             test_hmm.py (from hmm-m6a), test_pipeline.py, test_compare_glori.py
```
