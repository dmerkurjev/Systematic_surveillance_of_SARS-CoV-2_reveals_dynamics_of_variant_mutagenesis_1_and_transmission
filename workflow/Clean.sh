#!/bin/bash

# RNA-seq mini-pipeline: QC → alignment → counting

# -------------------- Setup --------------------


# Define project structure relative to current location
PROJECT_DIR="./data_pre_processing"
mkdir -p "${PROJECT_DIR}"/{raw,fastq,aligned,counts,logs,qc,bowtie2_index}


cd "${PROJECT_DIR}/raw"


# Group SRA run IDs by biological sample 
neg1=(SRX19322489)   # SRX19322489
neg2=(SRX19322492)   # SRX19322492 
pos1=(SRX19337747)    # SRX19337747
pos2=(SRX19347479)    # SRX7865902
pos3=(SRX7865901)    # SRX7865901
pos4=(SRX7865902)    # SRX7865902

# -------------------- Download & Convert --------------------

# Download .sra files
for r in "${neg1[@]}" "${neg2[@]}" "${pos1[@]}" "${pos2[@]}" "${pos1[@]}" "${pos2[@]}"; do
  prefetch "$r"
done

# Convert to gzipped FASTQ
for r in ${neg1[@]}" "${neg2[@]}" "${pos1[@]}" "${pos2[@]}" "${pos1[@]}" "${pos2[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done

# Concatenate per-sample FASTQs
cat "${neg1[@]/%/.fastq.gz}"  > ym.fastq.gz
cat "${neg2[@]/%/.fastq.gz}"  > yp.fastq.gz
cat "${pos1[@]/%/.fastq.gz}" > sm.fastq.gz
cat "${pos2[@]/%/.fastq.gz}" > sp.fastq.gz

# Move to fastq/ folder
mv n*.fastq.gz p*.fastq.gz ../fastq/

# -------------------- QC --------------------

cd ../fastq
fastqc ym.fastq.gz yp.gz sm.fastq.gz sp.fastq.gz \
  -o ../qc --threads 16

# -------------------- Alignment (hisat) --------------------

curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz

mkdir hisat2_index
cd hisat2_index
tar -xzf hg38.tar.gz

cd hisat2_index
SAMPLES=(ym,yp,sm,sp)

# Align each sample
for sample in "${SAMPLES[@]}"
do
  echo "Aligning ${sample}..."
  hisat2 -p 4 \
    -x hisat2_index/hg38/genome \
    -U fastq/${sample}.fastq.gz \
    2> logs/${sample}_hisat2.log | \
    samtools sort -@ 4 -o aligned/${sample}.bam
  samtools index aligned/${sample}.bam
  echo "${sample} alignment done."
done


# -------------------- Quantification (featureCounts) --------------------

cd ..
curl -L -o gencode.v48.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
gunzip -f gencode.v48.annotation.gtf.gz

featureCounts -T 16 -t exon -g gene_name \
  -a gencode.v48.annotation.gtf \
  -o counts/raw_counts_gene_sym.txt aligned/*.bam \
  &> logs/featureCounts_gene_sym.log

# Format counts matrix
{ printf "GeneSymbol\t"; head -n 2 counts/raw_counts_gene_sym.txt | tail -n 1 | cut -f7-; } > counts/final_counts_symbols.tsv
tail -n +3 counts/raw_counts_gene_sym.txt | \
  awk -v OFS="\t" '{ out=$1; for(i=7;i<=NF;i++) out=out OFS $i; print out }' >> Normalized_counts/final_counts_symbols.tsv

sed -i '' '1 s|aligned/||g; 1 s|\.bam||g' counts/final_counts_symbols.tsv

# Done
echo "Pipeline complete. Output saved in: ${PROJECT_DIR}/counts/final_counts_symbols.tsv"
