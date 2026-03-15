# QIIME2 Amplicon Pipeline

A pipeline for 16S rRNA amplicon sequencing analysis using QIIME2. Covers quality
control, denoising, taxonomic classification, phylogenetic tree construction,
diversity analysis, differential abundance analysis, and data export for downstream
analysis in R.

Developed based on course material by Dr. G. Rouskas.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Environment Setup](#2-environment-setup)
3. [Quality Control](#3-quality-control)
4. [Import Data into QIIME2](#4-import-data-into-qiime2)
5. [Adapter Trimming with Cutadapt](#5-adapter-trimming-with-cutadapt)
6. [Denoising with DADA2](#6-denoising-with-dada2)
7. [Taxonomic Classification with VSEARCH](#7-taxonomic-classification-with-vsearch)
8. [Phylogenetic Tree Construction](#8-phylogenetic-tree-construction)
9. [Diversity Analysis](#9-diversity-analysis)
10. [Rarefaction Analysis](#10-rarefaction-analysis)
11. [Differential Abundance Analysis with ANCOM](#11-differential-abundance-analysis-with-ancom)
12. [Export Data for R](#12-export-data-for-r)
13. [Acknowledgements](#13-acknowledgements)

---

## 1. Overview

QIIME2 (Quantitative Insights Into Microbial Ecology) is an open-source platform for
microbiome analysis. It works with its own file formats:

| Format | Description |
|--------|-------------|
| `.qza` | QIIME2 artifact — contains data and provenance information |
| `.qzv` | QIIME2 visualization — can be viewed at [view.qiime2.org](https://view.qiime2.org) |

All `.qzv` files produced in this pipeline can be interactively explored by dragging
them into [view.qiime2.org](https://view.qiime2.org).

---

## 2. Environment Setup

QIIME2 is installed and managed through conda. Load the required modules and activate
the environment before running any commands.

```bash
module load gcc/14.2.0 miniconda3/24.7.1
source $CONDA_PROFILE/conda.sh
conda activate qiime2-amplicon-2024.10
```

**Install QIIME2 (first time only):**

```bash
conda activate mamba_env

mamba env create -n qiime2-amplicon-2024.10 \
    --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml -y

conda activate qiime2-amplicon-2024.10
qiime info
```

---

## 3. Quality Control

FastQC generates a quality report per sample. MultiQC aggregates all FastQC reports
into a single summary report, making it easier to compare samples.

```bash
module load gcc/14.2.0 fastqc/0.12.1

# Run FastQC on all samples
fastqc -t 20 *.fastq.gz -o fastqc/
```

```bash
# Run MultiQC to summarise all FastQC reports
multiqc fastqc/ -o multiqc/
```

- `-t 20` uses 20 threads to speed up processing
- The FastQC HTML reports are saved in the `fastqc/` directory
- The MultiQC summary report is saved in the `multiqc/` directory

---

## 4. Import Data into QIIME2

Before any analysis, the raw FASTQ files must be imported into QIIME2 format (.qza).
This requires a manifest file, which is a tab-separated file that maps each sample
name to its forward and reverse FASTQ file paths.

**Manifest file format (`manifest.tsv`):**

```
sample-id    absolute-filepath                        direction
sample1      /path/to/sample1_R1.fastq.gz             forward
sample1      /path/to/sample1_R2.fastq.gz             reverse
sample2      /path/to/sample2_R1.fastq.gz             forward
sample2      /path/to/sample2_R2.fastq.gz             reverse
```

**Import the data:**

```bash
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest.tsv \
    --output-path demux.qza \
    --input-format PairedEndFastqManifestPhred33
```

**Visualise the imported data:**

```bash
qiime demux summarize \
    --i-data demux.qza \
    --o-visualization demux.qzv
```

Open `demux.qzv` at [view.qiime2.org](https://view.qiime2.org) to inspect read
quality and decide trimming parameters for the next step.

---

## 5. Adapter Trimming with Cutadapt

Cutadapt removes primer sequences from the reads. This step is important for amplicon
data because the primer sequences are not part of the biological sequence and must be
removed before denoising.

```bash
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux.qza \
    --p-cores 20 \
    --o-trimmed-sequences demux_trimmed.qza \
    --verbose
```

**Visualise the trimmed data:**

```bash
qiime demux summarize \
    --i-data demux_trimmed.qza \
    --o-visualization demux_trimmed.qzv
```

Compare `demux_trimmed.qzv` with `demux.qzv` to confirm that primers were removed
and to choose the truncation lengths for DADA2 in the next step.

---

## 6. Denoising with DADA2

DADA2 performs denoising of the paired-end reads. It corrects sequencing errors,
removes chimeric sequences, and produces a table of Amplicon Sequence Variants (ASVs).
ASVs are exact sequence variants, offering higher resolution than OTUs.

The truncation lengths (`--p-trunc-len-f` and `--p-trunc-len-r`) and trim lengths
(`--p-trim-left-f` and `--p-trim-left-r`) should be chosen based on the quality
plots from the previous step. Truncate where quality drops below Q25-Q30.

```bash
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux_trimmed.qza \
    --p-trim-left-f 17 \
    --p-trim-left-r 21 \
    --p-trunc-len-f 260 \
    --p-trunc-len-r 219 \
    --p-min-overlap 12 \
    --p-chimera-method consensus \
    --p-n-threads 20 \
    --o-table dada2_table.qza \
    --o-representative-sequences dada2_seqs.qza \
    --o-denoising-stats dada2_denoising_stats.qza
```

**DADA2 parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--p-trim-left-f` | 17 | Trim 17 bases from the start of forward reads (primer region) |
| `--p-trim-left-r` | 21 | Trim 21 bases from the start of reverse reads (primer region) |
| `--p-trunc-len-f` | 260 | Truncate forward reads at position 260 |
| `--p-trunc-len-r` | 219 | Truncate reverse reads at position 219 |
| `--p-min-overlap` | 12 | Minimum overlap between forward and reverse reads for merging |
| `--p-chimera-method` | consensus | Remove chimeric reads using the consensus method |

**Visualise denoising statistics:**

```bash
qiime metadata tabulate \
    --m-input-file dada2_denoising_stats.qza \
    --o-visualization dada2_denoising_stats.qzv
```

Inspect `dada2_denoising_stats.qzv` to check how many reads passed each filtering step.

---

## 7. Taxonomic Classification with VSEARCH

VSEARCH performs consensus-based taxonomic classification by comparing the ASV
sequences against a reference database (SILVA 132 at 97% identity).

```bash
qiime feature-classifier classify-consensus-vsearch \
    --i-query dada2_seqs.qza \
    --i-reference-reads silva-132-97-16s.qza \
    --i-reference-taxonomy silva-132-97-16s-taxonomy.qza \
    --p-maxaccepts 1 \
    --p-perc-identity 0.97 \
    --p-threads 20 \
    --o-classification dada2_vsearch_classification.qza \
    --o-search-results dada2_vsearch_hits.qza
```

**Parameter explanations:**

| Parameter | Description |
|-----------|-------------|
| `--p-maxaccepts 1` | Accept only the top hit for each query sequence |
| `--p-perc-identity 0.97` | Minimum 97% sequence identity to assign taxonomy |

**Visualise classification results:**

```bash
qiime metadata tabulate \
    --m-input-file dada2_vsearch_classification.qza \
    --o-visualization dada2_vsearch_classification.qzv
```

**Visualise taxonomy bar plots:**

```bash
qiime taxa barplot \
    --i-table dada2_table.qza \
    --i-taxonomy dada2_vsearch_classification.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization dada2_vsearch_classification_barplots.qzv
```

---

## 8. Phylogenetic Tree Construction

A phylogenetic tree is built from the ASV sequences using MAFFT for multiple sequence
alignment and FastTree for tree construction. The tree is required for phylogeny-based
diversity metrics such as Faith's PD and UniFrac distances.

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences dada2_seqs.qza \
    --output-dir phylogenetic_tree \
    --p-n-threads 10 \
    --verbose \
    &> phylogenetic_tree_generation.log
```

The output directory contains:

| File | Description |
|------|-------------|
| `aligned-dna-sequences.qza` | Multiple sequence alignment |
| `masked-aligned-dna-sequences.qza` | Masked alignment (removes uninformative positions) |
| `unrooted-tree.qza` | Unrooted phylogenetic tree |
| `rooted-tree.qza` | Midpoint-rooted tree — used in diversity analysis |

---

## 9. Diversity Analysis

Core diversity metrics are calculated using the rooted phylogenetic tree and the
feature table. The sampling depth normalizes all samples to the same number of reads
(rarefaction) and should be set to the minimum sequencing depth across all samples.

```bash
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny phylogenetic_tree/rooted-tree.qza \
    --i-table dada2_table.qza \
    --p-sampling-depth 21325 \
    --output-dir core-metrics \
    --m-metadata-file metadata.tsv \
    --p-n-jobs 10 \
    --verbose \
    &> core_metrics.log
```

**Metrics calculated:**

| Type | Metric | Description |
|------|--------|-------------|
| Alpha | Faith's PD | Phylogenetic diversity within a sample |
| Alpha | Shannon index | Species richness and evenness |
| Alpha | Observed features | Number of unique ASVs per sample |
| Beta | Weighted UniFrac | Phylogenetic distance between samples (abundance-weighted) |
| Beta | Unweighted UniFrac | Phylogenetic distance between samples (presence/absence) |
| Beta | Bray-Curtis | Compositional dissimilarity between samples |

---

## 10. Rarefaction Analysis

Alpha rarefaction curves show how species richness changes with sequencing depth.
They help determine whether the sequencing depth is sufficient to capture the full
diversity of each sample. A curve that plateaus indicates sufficient depth.

```bash
qiime diversity alpha-rarefaction \
    --i-table dada2_table.qza \
    --i-phylogeny phylogenetic_tree/rooted-tree.qza \
    --p-max-depth 21325 \
    --o-visualization rarefaction.qzv
```

- `--p-max-depth` should match the sampling depth used in the diversity analysis

---

## 11. Differential Abundance Analysis with ANCOM

ANCOM (Analysis of Composition of Microbiomes) identifies taxa that are differentially
abundant between groups. Because microbiome data is compositional, a pseudocount is
added before log-transformation to avoid issues with zero counts.

**Add pseudocount:**

```bash
qiime composition add-pseudocount \
    --i-table dada2_table.qza \
    --o-composition-table comp-feature-table.qza
```

**Run ANCOM:**

```bash
qiime composition ancom \
    --i-table comp-feature-table.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column group \
    --p-transform-function clr \
    --o-visualization ancom-group.qzv
```

- `--m-metadata-column group` specifies which column in the metadata file defines
  the groups to compare
- `--p-transform-function clr` applies centered log-ratio (CLR) transformation;
  `log` can be used as an alternative

---

## 12. Export Data for R

To perform downstream analysis in R (e.g. with phyloseq or vegan), the QIIME2
artifacts need to be exported to standard formats.

**Export the feature table (count table):**

```bash
qiime tools export \
    --input-path dada2_table.qza \
    --output-path export/table

# Convert the .biom file to .tsv
biom convert \
    -i export/table/feature-table.biom \
    -o export/table/feature-table.tsv \
    --to-tsv
```

**Export representative sequences:**

```bash
qiime tools export \
    --input-path dada2_seqs.qza \
    --output-path export/rep-seqs
```

**Export taxonomy:**

```bash
qiime tools export \
    --input-path dada2_vsearch_classification.qza \
    --output-path export/taxonomy
```

**Export phylogenetic tree:**

```bash
qiime tools export \
    --input-path phylogenetic_tree/rooted-tree.qza \
    --output-path export/tree
```

The exported files can be loaded into R using the phyloseq package for further
statistical analysis and visualisation.

---

## 13. Acknowledgements

This pipeline is based on course material from the **Metagenomics** module of the MSc program "Applied Bioinformatics" (AUTH/IHU).

Original scripts developed by **Dr. G. Rouskas**.
