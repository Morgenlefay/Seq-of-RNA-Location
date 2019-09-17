# APEX-Seq

## APEX-Seq
PMID: [31230715](https://www.ncbi.nlm.nih.gov/pubmed/31230715)

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
ls APEX/*fastq.gz |  while read id; do 
$fastqc -f fastq -o APEX/FastQC_1/./ ${id} &
done

multiqc APEX/FastQC_1 -o APEX/FastQC_1
```

### Step3-Hisat2

```bash
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367268_1.fastq.gz -2 APEX/SRR7367268_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nls-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367269_1.fastq.gz -2 APEX/SRR7367269_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nls-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367270_1.fastq.gz -2 APEX/SRR7367270_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nls-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367271_1.fastq.gz -2 APEX/SRR7367271_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nls-control-2.bam

hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367272_1.fastq.gz -2 APEX/SRR7367272_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nik-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367273_1.fastq.gz -2 APEX/SRR7367273_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nik-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367274_1.fastq.gz -2 APEX/SRR7367274_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nik-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367275_1.fastq.gz -2 APEX/SRR7367275_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nik-control-2.bam

hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367276_1.fastq.gz -2 APEX/SRR7367276_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Lma-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367277_1.fastq.gz -2 APEX/SRR7367277_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Lma-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367278_1.fastq.gz -2 APEX/SRR7367278_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Lma-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR8646402_1.fastq.gz -2 APEX/SRR8646402_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Lma-control-2.bam

hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367280_1.fastq.gz -2 APEX/SRR7367280_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nucpore-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367281_1.fastq.gz -2 APEX/SRR7367281_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nucpore-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367282_1.fastq.gz -2 APEX/SRR7367282_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nucpore-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367283_1.fastq.gz -2 APEX/SRR7367283_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nucpore-control-2.bam

hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367284_1.fastq.gz -2 APEX/SRR7367284_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nes-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367285_1.fastq.gz -2 APEX/SRR7367285_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nes-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367286_1.fastq.gz -2 APEX/SRR7367286_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nes-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367287_1.fastq.gz -2 APEX/SRR7367287_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nes-control-2.bam

hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367288_1.fastq.gz -2 APEX/SRR7367288_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Erm-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367289_1.fastq.gz -2 APEX/SRR7367289_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Erm-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367290_1.fastq.gz -2 APEX/SRR7367290_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Erm-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367291_1.fastq.gz -2 APEX/SRR7367291_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Erm-control-2.bam

hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367292_1.fastq.gz -2 APEX/SRR7367292_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Kdel-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367293_1.fastq.gz -2 APEX/SRR7367293_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Kdel-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367294_1.fastq.gz -2 APEX/SRR7367294_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Kdel-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367295_1.fastq.gz -2 APEX/SRR7367295_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Kdel-control-2.bam

hisat2 -t -p 8 -x Reference/index/mm10/genome -1 APEX/SRR7367296_1.fastq.gz -2 APEX/SRR7367296_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Omm-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367297_1.fastq.gz -2 APEX/SRR7367297_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Omm-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367298_1.fastq.gz -2 APEX/SRR7367298_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Omm-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367299_1.fastq.gz -2 APEX/SRR7367299_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Omm-control-2.bam

hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367300_1.fastq.gz -2 APEX/SRR7367300_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Mito-target-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367301_1.fastq.gz -2 APEX/SRR7367301_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Mito-target-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367302_1.fastq.gz -2 APEX/SRR7367302_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Mito-target-3.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367303_1.fastq.gz -2 APEX/SRR7367303_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Mito-control-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367304_1.fastq.gz -2 APEX/SRR7367304_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Mito-control-2.bam

hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367333_1.fastq.gz -2 APEX/SRR7367333_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nuclear-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367334_1.fastq.gz -2 APEX/SRR7367334_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Nuclear-2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367335_1.fastq.gz -2 APEX/SRR7367335_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Cytosol-1.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 APEX/SRR7367336_1.fastq.gz -2 APEX/SRR7367336_2.fastq.gz | samtools sort -@4 -O bam -o APEX/Cytosol-2.bam
```
### Step4-Stringtie

```bash
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nuclear_1/APEX_Nuclear_1.gtf APEX/Nuclear-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nuclear_2/APEX_Nuclear_2.gtf APEX/Nuclear-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Cytosol_1/APEX_Cytosol_1.gtf APEX/Cytosol-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Cytosol_2/APEX_Cytosol_2.gtf APEX/Cytosol-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nls_target_1/APEX_Nls_target_1.gtf APEX/Nls-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nls_target_2/APEX_Nls_target_2.gtf APEX/Nls-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nls_control_1/APEX_Nls_control_1.gtf APEX/Nls-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nls_control_2/APEX_Nls_control_2.gtf APEX/Nls-control-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nik_target_1/APEX_Nik_target_1.gtf APEX/Nik-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nik_target_2/APEX_Nik_target_2.gtf APEX/Nik-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nik_control_1/APEX_Nik_control_1.gtf APEX/Nik-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nik_control_2/APEX_Nik_control_2.gtf APEX/Nik-control-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Lma_target_1/APEX_Lma_target_1.gtf APEX/Lma-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Lma_target_2/APEX_Lma_target_2.gtf APEX/Lma-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Lma_control_1/APEX_Lma_control_1.gtf APEX/Lma-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Lma_control_2/APEX_Lma_control_2.gtf APEX/Lma-control-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nucpore_target_1/APEX_Nucpore_target_1.gtf APEX/Nucpore-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nucpore_target_2/APEX_Nucpore_target_2.gtf APEX/Nucpore-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nucpore_control_1/APEX_Nucpore_control_1.gtf APEX/Nucpore-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nucpore_control_2/APEX_Nucpore_control_2.gtf APEX/Nucpore-control-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nes_target_1/APEX_Nes_target_1.gtf APEX/Nes-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nes_target_2/APEX_Nes_target_2.gtf APEX/Nes-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nes_control_1/APEX_Nes_control_1.gtf APEX/Nes-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Nes_control_2/APEX_Nes_control_2.gtf APEX/Nes-control-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Erm_target_1/APEX_Erm_target_1.gtf APEX/Erm-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Erm_target_2/APEX_Erm_target_2.gtf APEX/Erm-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Erm_control_1/APEX_Erm_control_1.gtf APEX/Erm-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Erm_control_2/APEX_Erm_control_2.gtf APEX/Erm-control-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Kdel_target_1/APEX_Kdel_target_1.gtf APEX/Kdel-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Kdel_target_2/APEX_Kdel_target_2.gtf APEX/Kdel-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Kdel_control_1/APEX_Kdel_control_1.gtf APEX/Kdel-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Kdel_control_2/APEX_Kdel_control_2.gtf APEX/Kdel-control-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Omm_target_1/APEX_Omm_target_1.gtf APEX/Omm-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Omm_target_2/APEX_Omm_target_2.gtf APEX/Omm-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Omm_control_1/APEX_Omm_control_1.gtf APEX/Omm-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Omm_control_2/APEX_Omm_control_2.gtf APEX/Omm-control-2.bam

stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Mito_target_1/APEX_Mito_target_1.gtf APEX/Mito-target-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Mito_target_2/APEX_Mito_target_2.gtf APEX/Mito-target-2.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Mito_target_3/APEX_Mito_target_3.gtf APEX/Mito-target-3.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Mito_control_1/APEX_Mito_control_1.gtf APEX/Mito-control-1.bam
stringtie -e -B -p 8 -G Reference/GTF/gencode.v19.annotation.gtf -o ballgown/APEX_Mito_control_2/APEX_Mito_control_2.gtf APEX/Mito-control-2.bam
```
### Step5-Ballgown

```bash
###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(ballgown)
bg = ballgown(dataDir="ballgown", samplePattern = "APEX", meas='all')
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
### Step6-Bam2BigWig

```bash
chmod 775 bam2bigwig.sh
./bam2bigwig.sh APEX/Nuclear-1.bam
./bam2bigwig.sh APEX/Nuclear-2.bam
./bam2bigwig.sh APEX/Cytosol-1.bam
./bam2bigwig.sh APEX/Cytosol-2.bam
./bam2bigwig.sh APEX/Nls-target-1.bam
./bam2bigwig.sh APEX/Nls-target-2.bam
./bam2bigwig.sh APEX/Nls-control-1.bam
./bam2bigwig.sh APEX/Nls-control-2.bam
./bam2bigwig.sh APEX/Nik-target-1.bam
./bam2bigwig.sh APEX/Nik-target-2.bam
./bam2bigwig.sh APEX/Nik-control-1.bam
./bam2bigwig.sh APEX/Nik-control-2.bam
./bam2bigwig.sh APEX/Lma-target-1.bam
./bam2bigwig.sh APEX/Lma-target-2.bam
./bam2bigwig.sh APEX/Lma-control-1.bam
./bam2bigwig.sh APEX/Lma-control-2.bam
./bam2bigwig.sh APEX/Nucpore-target-1.bam
./bam2bigwig.sh APEX/Nucpore-target-2.bam
./bam2bigwig.sh APEX/Nucpore-control-1.bam
./bam2bigwig.sh APEX/Nucpore-control-2.bam
./bam2bigwig.sh APEX/Nes-target-1.bam
./bam2bigwig.sh APEX/Nes-target-2.bam
./bam2bigwig.sh APEX/Nes-control-1.bam
./bam2bigwig.sh APEX/Nes-control-2.bam
./bam2bigwig.sh APEX/Erm-target-1.bam
./bam2bigwig.sh APEX/Erm-target-2.bam
./bam2bigwig.sh APEX/Erm-control-1.bam
./bam2bigwig.sh APEX/Erm-control-2.bam
./bam2bigwig.sh APEX/Kdel-target-1.bam
./bam2bigwig.sh APEX/Kdel-target-2.bam
./bam2bigwig.sh APEX/Kdel-control-1.bam
./bam2bigwig.sh APEX/Kdel-control-2.bam
./bam2bigwig.sh APEX/Omm-target-1.bam
./bam2bigwig.sh APEX/Omm-target-2.bam
./bam2bigwig.sh APEX/Omm-control-1.bam
./bam2bigwig.sh APEX/Omm-control-2.bam
./bam2bigwig.sh APEX/Mito-target-1.bam
./bam2bigwig.sh APEX/Mito-target-2.bam
./bam2bigwig.sh APEX/Mito-target-3.bam
./bam2bigwig.sh APEX/Mito-control-1.bam
./bam2bigwig.sh APEX/Mito-control-2.bam
```
