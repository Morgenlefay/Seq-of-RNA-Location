## cheRNA
### Data
PMID: [26257179](https://www.ncbi.nlm.nih.gov/pubmed/26257179)

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
ls cheRNA/*fastq.gz |  while read id; do 
$fastqc -f fastq -o cheRNA/FastQC_1/./ ${id} &
done
multiqc cheRNA/FastQC_1 -o cheRNA/FastQC_1
```
### Step3-Hisat2
```bash
hisat2 -t -p 8 -x Reference/index/hg19/genome -U cheRNA/SRR1824492.fastq.gz | samtools sort -@4 -O bam -o cheRNA/SNE_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -U cheRNA/SRR1824493.fastq.gz | samtools sort -@4 -O bam -o cheRNA/SNE_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -U cheRNA/SRR1824494.fastq.gz | samtools sort -@4 -O bam -o cheRNA/SNE_3.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -U cheRNA/SRR1824495.fastq.gz | samtools sort -@4 -O bam -o cheRNA/CPE_1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -U cheRNA/SRR1824496.fastq.gz | samtools sort -@4 -O bam -o cheRNA/CPE_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -U cheRNA/SRR1824497.fastq.gz | samtools sort -@4 -O bam -o cheRNA/CPE_3.bam
```
### Step4-Bam2BigWig
```bash
chmod 775 bam2bigwig.sh
./bam2bigwig.sh cheRNA/SNE_1.bam
./bam2bigwig.sh cheRNA/SNE_2.bam
./bam2bigwig.sh cheRNA/SNE_3.bam
./bam2bigwig.sh cheRNA/CPE_1.bam
./bam2bigwig.sh cheRNA/CPE_2.bam
./bam2bigwig.sh cheRNA/CPE_3.bam
```
### Step5-Stringtie
```bash
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/cheRNA_SNE_1/cheRNA_SNE_1.gtf cheRNA/SNE_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/cheRNA_SNE_2/cheRNA_SNE_2.gtf cheRNA/SNE_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/cheRNA_SNE_3/cheRNA_SNE_3.gtf cheRNA/SNE_3.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/cheRNA_CPE_1/cheRNA_CPE_1.gtf cheRNA/CPE_1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/cheRNA_CPE_2/cheRNA_CPE_2.gtf cheRNA/CPE_2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/cheRNA_CPE_3/cheRNA_CPE_3.gtf cheRNA/CPE_3.bam
```

### Step6-Ballgown
```R
###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(ballgown)
bg = ballgown(dataDir="ballgown", samplePattern = "cheRNA", meas='all')
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
### Step7-Heatmap
TBtools (https://github.com/CJ-Chen/TBtools)






