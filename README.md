# ChIP-seq/DAP-seq analysis code and data

This repository serves as a collection of publicly available ChIP-seq/DAP-seq datasets (moslty in maize), provides codes / pipeline to systematically analyze them as well as QC and summarized data output.

Briefly:
- reads were trimmed by Trim Galore! and Cutadapt and aligned to the B73 maize reference genome using BWA.
- Picard was used to merge alignments from multiple libraries of the same sample and then to mark PCR duplicates
- Further filtering was employed to remove reads marked as:
  - duplicates
  - non-primary alignments
  - unmapped
  - mapped to multiple locations
  - having an insert size > 2kb
  - mapped with abnormal paired-end signatures:
    - only one read out of a pair mapped
    - two reads mapped to two different chromosomes
    - two reads mapped in abnormal orientation
- MACS2 was used to call broad and narrow peaks ("--gsize=2.1e9 --broad-cutoff=0.1") and calculate the FRiP scores.
- HOMER was used to annotate peaks relative to gene features.
- Genome-wide IP enrichment relative to control was calculated using deepTools
- Strand cross-correlation peak and ChIP-seq quality measures including NSC and RSC were calculated using phantompeakqualtools.

See [this table](/data/01.cfg.tsv) for a list of collected datasets.

[Check here](output.md) for a walk through of output files.