#!/bin/bash

# RNA-seq mini-pipeline: QC → alignment → counting

# -------------------- Setup --------------------


# Define project structure relative to current location
PROJECT_DIR="./data_pre_processing"
mkdir -p "${PROJECT_DIR}"/{raw,fastq,aligned,counts,logs,qc,bowtie2_index}


cd "${PROJECT_DIR}/raw"


# Group SRA run IDs by biological sample 
norm1=(SRX25770698)   # SRX25770698
norm2=(SRX25770699)   # SRX25770699
norm3=(SRX25770700)    # SRX25770700
tum1=(SRX25770708)    # SRX25770708
tum2=(SRX25770709)    # SRX25770709
tum3=(SRX25770710)    # SRX25770710

# -------------------- Download & Convert --------------------

# Download .sra files
for r in "${norm1[@]}" "${norm2[@]}" "${norm3[@]}" "${tum1[@]}" "${tum2[@]}" "${tum3[@]}"; do
  prefetch "$r"
done

# Convert to gzipped FASTQ
for r in "${norm1[@]}" "${norm2[@]}" "${norm3[@]}" "${tum1[@]}" "${tum2[@]}" "${tum3[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done

# Concatenate per-sample FASTQs
cat "${norm1[@]/%/.fastq.gz}"  > norm1.fastq.gz
cat "${norm2[@]/%/.fastq.gz}"  > norm2.fastq.gz
cat "${norm3[@]/%/.fastq.gz}" > norm3.fastq.gz
cat "${tum1[@]/%/.fastq.gz}" > tum1.fastq.gz
cat "${tum2[@]/%/.fastq.gz}" > tum2.fastq.gz
cat "${tum3[@]/%/.fastq.gz}" > tum3.fastq.gz

# Move to fastq/ folder
mv n*.fastq.gz p*.fastq.gz ../fastq/

# -------------------- QC --------------------

cd ../fastq
fastqc norm1.fastq.gz norm2.gz norm3.fastq.gz tum1.fastq.gz tum2.fastq.gz tum3.fastq.gz \
  -o ../qc --threads 16

# -------------------- Alignment (hisat) --------------------

curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/bowtie2/data/hg38.tar.gz

mkdir bowtie2_index
cd bowtie2_index
tar -xzf hg38.tar.gz

cd hisat2_index
SAMPLES=(pos1,pos2,neg1,neg2,neg3,neg4)

# Align each sample
for sample in "${SAMPLES[@]}"
do
  echo "Aligning ${sample}..."
  hisat2 -p 4 \
    -x bowtie2_index/hg38/genome \
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
