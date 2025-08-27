 SARS-CoV-2 RNA-seq — Preprocessing & Analysis

This repository contains an end-to-end workflow for RNA-seq analysis of early passage and senescent BJ fibroblasts in the absence or continued presence of doxycycline (20nm).
It includes two main components:

Preprocessing (Clean.sh)
Downloads SRA runs, builds 4 sample FASTQs, performs QC, alignment, and outputs gene transcription level.

Analysis (Differential_expression.sh)
Performs differential expression analysis using edgeR.

Data

SRA BioProject: 39738001 
RNAseq and sequence variant analysis for 63036 human nasopharyngeal swabs for SARS-CoV-2, with 1928 negative swabs and 2886 synthetic controls.
(part2) Additional RNAseq and sequence variant analysis for 4 human nasopharyngeal swabs for SARS-CoV-2, with 2 negative swabs.

Sample 1: GSM7036174 (SRX19322489)
tissue: human nasopharyngeal swab
sample type: negative sample
pipeline.version: preclinical

Sample 2: GSM4395597 (SRX19322492)
tissue: human nasopharyngeal swab
sample type: negative sample
pipeline.version: preclinical

Sample 3: GSM7039241 (SRX19337747)
tissue: human nasopharyngeal swab
sample type: positive sample
pipeline.version: preclinical

Sample 4: GSM7039810 (SRX19347479)
tissue: human nasopharyngeal swab
sample type: positive sample
pipeline.version: preclinical

Sample 5: GSM7039241 (GSM7257036)
tissue: human nasopharyngeal swab
sample type: positive sample
pipeline.version: preclinical

Sample 6: GSM7039810 (GSM7257037)
tissue: human nasopharyngeal swab
sample type: positive sample
pipeline.version: preclinical

Reference files:

Transcript annotation: clinical_V2_sparsq_V1p2_amplicons.gtf
Alignment index: hg38
Related publication:
PMID: 39738001

Directory Output

The preprocessing script creates a structured working directory data_pre_processing/ containing:

raw/      # downloaded and processed SRA FASTQs  
fastq/    # concatenated FASTQs (one per sample)  
aligned/  # hisat-aligned BAMs  
counts/   # gene-level count matrix   
qc/       # FastQC output  
bowtie2_index/  # reference index for alignment  


You can skip raw FASTQ processing and go straight to analysis using data_RNAcounts/final_counts_symbols.tsv.
