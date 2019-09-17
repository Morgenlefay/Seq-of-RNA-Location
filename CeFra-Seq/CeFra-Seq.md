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

fastqc=/home/morgenlefay/miniconda2/pkgs/fastqc-0.11.8-1/bin/fastqc
ls CeFra_seq/SK-N-SH/*fastq.gz |  while read id; do 
$fastqc -f fastq -o CeFra_seq/FastQC_1/./ ${id} &
done
multiqc CeFra_seq/FastQC_1 -o CeFra_seq/FastQC_1
```

### Step3-Hisat2

```bash
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HepG2/Cytosol_1_1.fastq.gz -2 CeFra_seq/HepG2/Cytosol_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HepG2/Cytosol_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HepG2/Cytosol_2_1.fastq.gz -2 CeFra_seq/HepG2/Cytosol_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HepG2/Cytosol_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HepG2/Insoluble_1_1.fastq.gz -2 CeFra_seq/HepG2/Insoluble_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HepG2/Insoluble_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HepG2/Insoluble_2_1.fastq.gz -2 CeFra_seq/HepG2/Insoluble_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HepG2/Insoluble_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HepG2/Membrane_1_1.fastq.gz -2 CeFra_seq/HepG2/Membrane_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HepG2/Membrane_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HepG2/Membrane_2_1.fastq.gz -2 CeFra_seq/HepG2/Membrane_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HepG2/Membrane_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HepG2/Nuclear_1_1.fastq.gz -2 CeFra_seq/HepG2/Nuclear_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HepG2/Nuclear_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HepG2/Nuclear_2_1.fastq.gz -2 CeFra_seq/HepG2/Nuclear_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HepG2/Nuclear_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/K562/Cytosol_1_1.fastq.gz -2 CeFra_seq/K562/Cytosol_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/K562/Cytosol_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/K562/Cytosol_2_1.fastq.gz -2 CeFra_seq/K562/Cytosol_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/K562/Cytosol_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/K562/Insoluble_1_1.fastq.gz -2 CeFra_seq/K562/Insoluble_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/K562/Insoluble_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/K562/Insoluble_2_1.fastq.gz -2 CeFra_seq/K562/Insoluble_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/K562/Insoluble_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/K562/Membrane_1_1.fastq.gz -2 CeFra_seq/K562/Membrane_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/K562/Membrane_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/K562/Membrane_2_1.fastq.gz -2 CeFra_seq/K562/Membrane_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/K562/Membrane_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/K562/Nuclear_1_1.fastq.gz -2 CeFra_seq/K562/Nuclear_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/K562/Nuclear_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/K562/Nuclear_2_1.fastq.gz -2 CeFra_seq/K562/Nuclear_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/K562/Nuclear_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HeLa-S3/Cytosol_1_1.fastq.gz -2 CeFra_seq/HeLa-S3/Cytosol_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HeLa-S3/Cytosol_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HeLa-S3/Cytosol_2_1.fastq.gz -2 CeFra_seq/HeLa-S3/Cytosol_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HeLa-S3/Cytosol_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HeLa-S3/Nuclear_1_1.fastq.gz -2 CeFra_seq/HeLa-S3/Nuclear_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HeLa-S3/Nuclear_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/HeLa-S3/Nuclear_2_1.fastq.gz -2 CeFra_seq/HeLa-S3/Nuclear_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/HeLa-S3/Nuclear_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/SK-N-SH/Cytosol_1_1.fastq.gz -2 CeFra_seq/SK-N-SH/Cytosol_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/SK-N-SH/Cytosol_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/SK-N-SH/Cytosol_2_1.fastq.gz -2 CeFra_seq/SK-N-SH/Cytosol_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/SK-N-SH/Cytosol_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/SK-N-SH/Nuclear_1_1.fastq.gz -2 CeFra_seq/SK-N-SH/Nuclear_1_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/SK-N-SH/Nuclear_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 CeFra_seq/SK-N-SH/Nuclear_2_1.fastq.gz -2 CeFra_seq/SK-N-SH/Nuclear_2_2.fastq.gz | samtools sort -@4 -O bam -o CeFra_seq/SK-N-SH/Nuclear_2.bam
```

### Step4-Stringtie

```bash
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HepG2_Cytosol_1/CeFra_HepG2_Cytosol_1.gtf CeFra_seq/HepG2/Cytosol_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HepG2_Cytosol_2/CeFra_HepG2_Cytosol_2.gtf CeFra_seq/HepG2/Cytosol_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HepG2_Insoluble_1/CeFra_HepG2_Insoluble_1.gtf CeFra_seq/HepG2/Insoluble_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HepG2_Insoluble_2/CeFra_HepG2_Insoluble_2.gtf CeFra_seq/HepG2/Insoluble_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HepG2_Membrane_1/CeFra_HepG2_Membrane_1.gtf CeFra_seq/HepG2/Membrane_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HepG2_Membrane_2/CeFra_HepG2_Membrane_2.gtf CeFra_seq/HepG2/Membrane_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HepG2_Nuclear_1/CeFra_HepG2_Nuclear_1.gtf CeFra_seq/HepG2/Nuclear_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HepG2_Nuclear_2/CeFra_HepG2_Nuclear_2.gtf CeFra_seq/HepG2/Nuclear_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_K562_Cytosol_1/CeFra_K562_Cytosol_1.gtf CeFra_seq/K562/Cytosol_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_K562_Cytosol_2/CeFra_K562_Cytosol_2.gtf CeFra_seq/K562/Cytosol_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_K562_Insoluble_1/CeFra_K562_Insoluble_1.gtf CeFra_seq/K562/Insoluble_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_K562_Insoluble_2/CeFra_K562_Insoluble_2.gtf CeFra_seq/K562/Insoluble_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_K562_Membrane_1/CeFra_K562_Membrane_1.gtf CeFra_seq/K562/Membrane_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_K562_Membrane_2/CeFra_K562_Membrane_2.gtf CeFra_seq/K562/Membrane_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_K562_Nuclear_1/CeFra_K562_Nuclear_1.gtf CeFra_seq/K562/Nuclear_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_K562_Nuclear_2/CeFra_K562_Nuclear_2.gtf CeFra_seq/K562/Nuclear_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HeLa-S3_Cytosol_1/CeFra_HeLa-S3_Cytosol_1.gtf CeFra_seq/HeLa-S3/Cytosol_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HeLa-S3_Cytosol_2/CeFra_HeLa-S3_Cytosol_2.gtf CeFra_seq/HeLa-S3/Cytosol_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HeLa-S3_Nuclear_1/CeFra_HeLa-S3_Nuclear_1.gtf CeFra_seq/HeLa-S3/Nuclear_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_HeLa-S3_Nuclear_2/CeFra_HeLa-S3_Nuclear_2.gtf CeFra_seq/HeLa-S3/Nuclear_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_SK-N-SH_Cytosol_1/CeFra_SK-N-SH_Cytosol_1.gtf CeFra_seq/SK-N-SH/Cytosol_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_SK-N-SH_Cytosol_2/CeFra_SK-N-SH_Cytosol_2.gtf CeFra_seq/SK-N-SH/Cytosol_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_SK-N-SH_Nuclear_1/CeFra_SK-N-SH_Nuclear_1.gtf CeFra_seq/SK-N-SH/Nuclear_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/CeFra_SK-N-SH_Nuclear_2/CeFra_SK-N-SH_Nuclear_2.gtf CeFra_seq/SK-N-SH/Nuclear_2.bam
```
### Step5-Ballgown

```R
###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(ballgown)
bg = ballgown(dataDir="ballgown", samplePattern = "CeFra", meas='all')
gene_expression = gexpr(bg)
write.csv(gene_expression, file="gene_expression.csv",row.names = T)

###FPKM list
Gene_expression <- read.csv("gene_expression.csv")
ENSEMBL <- gsub("\\.\\d*", "", Gene_expression$Geneid)
write.csv(ENSEMBL,file="ENSEMBL.csv")
Gene_expression<- read.csv("gene_expression.csv",header = TRUE)
index<-duplicated(Gene_expression$Geneid)
Gene_expression<-Gene_expression[!index,]
write.csv(Gene_expression, file="gene_expression.csv",row.names = F)

###Genename list
library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, Geneid = ensembl_gene_id, ID = external_gene_name)

###Merge
Gene_expression <- read.csv("gene_expression.csv")
Gene_expression <- merge(Gene_expression,t2g)
write.csv(Gene_expression, file="Gene_expression.csv",row.names = F)
```
