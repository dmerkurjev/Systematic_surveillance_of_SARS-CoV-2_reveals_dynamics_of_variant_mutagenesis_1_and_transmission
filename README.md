CRISP-based mouse model of Vhl-deficient clear cell kidney cancer — Preprocessing & Analysis

This repository contains an end-to-end workflow for RNA-Seq and CRISP analysis of 3 normal and 3 tumor kidney samples.
It includes two main components:

Preprocessing (Clean.sh)
Downloads SRA runs, builds 4 sample FASTQs, performs QC, alignment, and outputs gene expression level.

Analysis (CRISP1.R and CRISP2.R)
Performs RNA-Seq and CRISP analysis using edgeR.

Data

SRA BioProject: PRJNA1150068
RNAseq and CRISP analysis for 3 normal and 3 tumor kidney samples.

Sample 1: GSM8473211 (SRX25770698)
tissue: kidney
group: Normal kidney

Sample 2: GSM8473212 (SRX25770699)
tissue: kidney
group: Normal kidney

Sample 3: GSM8473213 (SRX25770700)
tissue: kidney
group: Normal kidney

Sample 4: GSM7275379 (SRX25770708)
tissue: kidney
group: Tumor

Sample 5: GSM7275380 (SRX25770709)
tissue: kidney
group: Tumor

Sample 6: GSM7275381 (SRX25770710)
tissue: kidney
group: Tumor

Reference files:

Transcript annotation: tissue: gencode.vM38.annotation.gtf
Alignment index: mm10
Related publication:
PMID: 39365820

Directory Output

The preprocessing script creates a structured working directory data_pre_processing/ containing:

raw/      # downloaded and processed SRA FASTQs  
fastq/    # concatenated FASTQs (one per sample)  
counts/   # gene-level count matrix   
qc/       # FastQC output  
bowtie2_index/  # reference index for alignment  


You can skip raw FASTQ processing and go straight to analysis using data_RNAcounts/final_counts_symbols.tsv.
