###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(ballgown)
bg = ballgown(dataDir="ballgown", samplePattern = "APEX_Nls", meas='all')
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

###Dubr
Gene_expression <- read.csv("Gene_expression.csv")
Dubr <- read.csv("Dubr.csv")
Dubr <- merge(Dubr,Gene_expression)
write.csv(Dubr, file="Dubr_expression.csv",row.names = F)


