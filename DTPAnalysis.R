install.packages("DataCombine")
library(DataCombine)
library(readxl)
library(dplyr)
library(readr)
library(wPerm)
library(tidyr)
library(data.table)
DTP_NCI60_ZSCORE <- read_excel("DTP_NCI60_ZSCORE.xlsx")
geneExpression<-read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

Camptothecin<-DTP_NCI60_ZSCORE[DTP_NCI60_ZSCORE$`https://discover.nci.nih.gov/cellminer/`==c("Camptothecin", "Drug name"),]
Camptothecin<-Camptothecin[-c(1, 2,3),]
Camptothecin<-Camptothecin[,-c(3:6)]
Camp<-data.frame(lapply(Camptothecin, function(x){
    sub(".*:","", x)
}))

CampMerge<-as.data.frame(t(Camp))

rownames(CampMerge)<-c()
CampMerge<-CampMerge[-c(1,2,63,64),]

Genes<-as.list(unique(geneExpression$GENE_NAME))
