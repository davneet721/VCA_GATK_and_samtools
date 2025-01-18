# Script to call germline variants in a human WGS paired end reads 2 X 100bp
Following GATK4 best practices workflow: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

## Data files used:

```bash
# download FASTQ files
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

# download reference files
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# download known sites files for BQSR from GATK resource bundle
# wget -P /Users/davneetkaur/Desktop/Personal_Git_Projects/VCA_projects/Supporting_Files/hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
# wget -P /Users/davneetkaur/Desktop/Personal_Git_Projects/VCA_projects/Supporting_Files/hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
```
# Summary of further steps in the shell file:
```bash
###################################################### VARIANT CALLING STEPS ####################################################################

# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------

# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------

# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------
