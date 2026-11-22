"""hmm_m6a_fastq: m6A detection from clean GLORI-seq FASTQ (GLORI-tools references) with a per-nucleotide HMM."""
from hmm_m6a_fastq.hmm import (
    HMMParams,
    forward_backward,
    posterior_states,
    viterbi,
    emission_logprob,
)
from hmm_m6a_fastq.train import baum_welch
from hmm_m6a_fastq.pileup import pileup_bam, Region, AdenosineRecord
from hmm_m6a_fastq.calling import call_sites
from hmm_m6a_fastq.align import GloriRefs, align_fastq
from hmm_m6a_fastq.annotation import annotate_sites, load_tbl2

__version__ = "0.1.0"
__all__ = [
    "HMMParams",
    "forward_backward",
    "posterior_states",
    "viterbi",
    "emission_logprob",
    "baum_welch",
    "pileup_bam",
    "Region",
    "AdenosineRecord",
    "call_sites",
    "align_fastq",
    "GloriRefs",
    "annotate_sites",
    "load_tbl2",
    "__version__",
]
