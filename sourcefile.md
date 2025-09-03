# To download the data:

# Create environment
conda create -n Systematic_surveillance_of_SARS-CoV-2 -c bioconda -c conda-forge \
  sra-tools fastqc multiqc hisat2 samtools trimmomatic subread -y
conda activate Systematic_surveillance_of_SARS-CoV-2

# Folder setup
mkdir -p ~/0_Systematic_surveillance_of_SARS-CoV-2/{data,fastq,trimmed,aligned,counts,logs,qc}
cd ~/0_Systematic_surveillance_of_SARS-CoV-2/data

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR011/11255/SRR11255824/SRR11255824.lite.1 # negative sample nasal swab SRX20174770

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR011/11255/SRR11255825/SRR11255825.lite.1 # negative sample nasal swab SRX20174771

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR024/24387/SRR24387765/SRR24387765.lite.1 # negative sample nasal swab SRX20176102

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR024/24387/SRR24387929/SRR24387929.lite.1 # negative sample nasal swab SRX20176103

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos6/sra-pub-zq-40/SRR024/24387/SRR24387928/SRR24387928.lite.1 # positive sample nasal swab SRX20176104

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR011/11255/SRR11255827/SRR11255827.lite.1 # positive sample nasal swab SRX20176105

for r in "${SRR[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done
fastq-dump --split-files *lite.1 # convert to fastq
