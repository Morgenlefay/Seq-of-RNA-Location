## CeFra-Seq
### Methods and Data
PMID: [28579403](https://www.ncbi.nlm.nih.gov/pubmed/28579403)

PMID: [29079635](https://www.ncbi.nlm.nih.gov/pubmed/29079635)
### Step1-Set conda environment

```bash
conda create -n QC
conda activate QC
conda install -c bioconda fastqc
conda install -c bioconda rseqc
conda install -c bioconda trim-galore
conda install -c bioconda multiqc 
conda deactivate

conda create -n RNAseq
conda activate RNAseq
conda install -c bioconda hisat2
conda install -c bioconda bowtie2
conda install -c bioconda samtools
conda install -c bioconda bedtools 
conda install -c bioconda subread 
conda install -c bioconda htseq
conda install -c bioconda stringtie
conda install -c bioconda gffcompare
conda install -c bioconda rsem
conda install -c bioconda cufflinks
conda install -c bioconda tophat
conda deactivate
```
### Step2-Preprocess

#### FastQC and MultiQC of Raw Data

```bash
fastqc=/home/morgenlefay/miniconda2/pkgs/fastqc-0.11.8-1/bin/fastqc
ls CeFra_seq/HepG2/*fastq.gz |  while read id; do 
$fastqc -f fastq -o CeFra_seq/FastQC_1/./ ${id} &
done
multiqc CeFra_seq/FastQC_1 -o CeFra_seq/FastQC_1

fastqc=/home/morgenlefay/miniconda2/pkgs/fastqc-0.11.8-1/bin/fastqc
ls CeFra_seq/HeLa-S3/*fastq.gz |  while read id; do 
$fastqc -f fastq -o CeFra_seq/FastQC_1/./ ${id} &
done
multiqc CeFra_seq/FastQC_1 -o CeFra_seq/FastQC_1

fastqc=/home/morgenlefay/miniconda2/pkgs/fastqc-0.11.8-1/bin/fastqc
ls CeFra_seq/K562/*fastq.gz |  while read id; do 
$fastqc -f fastq -o CeFra_seq/FastQC_1/./ ${id} &
done
multiqc CeFra_seq/FastQC_1 -o CeFra_seq/FastQC_1
```

