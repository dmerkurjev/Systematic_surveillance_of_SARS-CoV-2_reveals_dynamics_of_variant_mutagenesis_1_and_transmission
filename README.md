SARS-CoV-2 RNA-seq — Preprocessing & Analysis

This repository contains an end-to-end workflow for RNA-seq analysis of early passage and senescent BJ fibroblasts in the absence or continued presence of doxycycline (20nm).
It includes two main components:

Preprocessing (Clean.sh)
Downloads SRA runs, builds 4 sample FASTQs, performs QC, alignment, and outputs gene transcription level.

Analysis (Differential_expression.sh)
Performs differential expression analysis using edgeR.

Data

SRA BioProject: 39738001 
RNAseq and sequence variant analysis for 6 human nasopharyngeal swabs for SARS-CoV-2, with 2 positive swabs and 4 negative swabs.

Sample 1: GSM7275266 (SRX20174770)
tissue: human nasopharyngeal swab
sample type: negative sample
pipeline.version: preclinical

Sample 2: GSM7275267 (SRX20174771)
tissue: human nasopharyngeal swab
sample type: negative sample
pipeline.version: preclinical

Sample 3: GSM7275378 (SRX20176102)
tissue: human nasopharyngeal swab
sample type: positive sample
pipeline.version: preclinical

Sample 4: GSM7275379 (SRX20176103)
tissue: human nasopharyngeal swab
sample type: positive sample
pipeline.version: preclinical

Sample 5: GSM7275380 (SRX20176104)
tissue: human nasopharyngeal swab
sample type: positive sample
pipeline.version: preclinical

Sample 6: GSM7275381 (SRX20176105)
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
aligned/  # bowtie2-aligned BAMs  
counts/   # gene-level count matrix   
qc/       # FastQC output  
bowtie2_index/  # reference index for alignment  


You can skip raw FASTQ processing and go straight to analysis using data_RNAcounts/final_counts_symbols.tsv.
