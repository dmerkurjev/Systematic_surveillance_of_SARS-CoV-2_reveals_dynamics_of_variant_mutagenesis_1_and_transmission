# To download the data:

# Create environment
conda create -n Toward_a_CRISPR-based_mouse_model_of_Vhl-deficient_clear_cell_kidney_cancer -c bioconda -c conda-forge \
  sra-tools fastqc multiqc hisat2 samtools trimmomatic subread -y
conda activate Toward_a_CRISPR-based_mouse_model_of_Vhl-deficient_clear_cell_kidney_cancer

# Folder setup
mkdir -p ~/0_Toward_a_CRISPR-based_mouse_model_of_Vhl-deficient_clear_cell_kidney_cancer/{data,fastq,trimmed,aligned,counts,logs,qc}
cd ~/0_Toward_a_CRISPR-based_mouse_model_of_Vhl-deficient_clear_cell_kidney_cancer/data

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR011/11255/SRR11255824/SRR11255824.lite.1 # normal kidney group 1 SRX25770699

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR011/11255/SRR11255825/SRR11255825.lite.1 # normal kidney group 2 SRX25770700

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR024/24387/SRR24387765/SRR24387765.lite.1 # normal kidney group 3 SRX25770701

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR024/24387/SRR24387929/SRR24387929.lite.1 # tumor kidney group 1 SRX25770708

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR024/24387/SRR24387928/SRR24387928.lite.1 # tumor kidney group 2 SRX25770709

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR011/11255/SRR11255827/SRR11255827.lite.1 # tumor kidney group 3 SRX25770710

for r in "${SRR[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done
fastq-dump --split-files *lite.1 # convert to fastq
