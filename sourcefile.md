# To download the data:

# Create environment
conda create -n Toward_a_CRISPR-based_mouse_model_of_Vhl-deficient_clear_cell_kidney_cancer -c bioconda -c conda-forge \
  sra-tools fastqc multiqc hisat2 samtools trimmomatic subread -y
conda activate Toward_a_CRISPR-based_mouse_model_of_Vhl-deficient_clear_cell_kidney_cancer

# Folder setup
mkdir -p ~/0_Toward_a_CRISPR-based_mouse_model_of_Vhl-deficient_clear_cell_kidney_cancer/{data,fastq,trimmed,aligned,counts,logs,qc}
cd ~/0_Toward_a_CRISPR-based_mouse_model_of_Vhl-deficient_clear_cell_kidney_cancer/data

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR030/30310/SRR30310278/SRR30310278.lite.1 # normal kidney group 1 SRX25770699

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR030/30310/SRR30310277/SRR30310277.lite.1 # normal kidney group 2 SRX25770700

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR030/30310/SRR30310276/SRR30310276.lite.1 # normal kidney group 3 SRX25770701

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR030/30310/SRR30310268/SRR30310268.lite.1 # tumor kidney group 1 SRX25770708

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR030/30310/SRR30310267/SRR30310267.lite.1 # tumor kidney group 2 SRX25770709

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR030/30310/SRR30310266/SRR30310266.lite.1 # tumor kidney group 3 SRX25770710

for r in "${SRR[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done
fastq-dump --split-files *lite.1 # convert to fastq
