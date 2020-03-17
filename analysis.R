library(dplyr)
library(readr)
library(wPerm)
library(data.table)
library(parallel)


cellLines<-readxl::read_xlsx("GDSC2_fitted_dose_response_25Feb20.xlsx")
cellLines<-cellLines%>%
    select("CELL_LINE_NAME", "DRUG_NAME", "LN_IC50")

geneExpression<-read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

geneExpression<-geneExpression%>%
    select(-"SAMPLE_ID", -"REGULATION")


getAllDrugs<-function(){
    cellLines<-readxl::read_xlsx("GDSC2_fitted_dose_response_25Feb20.xlsx")
    cellLines<-cellLines%>%
        select("CELL_LINE_NAME", "DRUG_NAME", "LN_IC50")
    geneExpression<-read_tsv("CosmicCLP_CompleteGeneExpression.tsv")
    geneExpression<-geneExpression%>%
        select(-"SAMPLE_ID", -"REGULATION")
    lapply(drugs, drugCorTest, geneFile=geneExpression, cellFile=cellLines)
}

drugCorTest<-function(geneFile, cellFile, drug){
    name<-drug
    selectLines<-cellLines %>%
        filter(DRUG_NAME == name)
    raw<-lapply(genesList, getCorr, geneFile=geneExpression, cellFile=cellLines)
    cleanup(name, raw)
}

getCorr<-function(geneFile, cellFile, geneName){
    selectGenes<-geneExpression %>%
        filter(GENE_NAME== geneName)
    
    merged<-left_join(selectLines,selectGenes, by=c("CELL_LINE_NAME"="SAMPLE_NAME"))
    merged<-na.omit(merged)
    x<-merged$LN_IC50
    y<-merged$Z_SCORE
    cor.test(x,y, method="spearman")
}

getPValue<- function(geneFile, cellFile, geneName){
    selectLines<-cellLines %>%
        filter(DRUG_NAME == "Cytarabine")
    selectGenes<-geneExpression %>%
        filter(GENE_NAME== geneName)
    
    merged<-left_join(selectLines,selectGenes, by=c("CELL_LINE_NAME"="SAMPLE_NAME"))
    merged<-na.omit(merged)
    x<-merged$LN_IC50
    y<-merged$Z_SCORE
    cor.test(x,y, method="spearman")
}


genes<-unique(geneExpression$GENE_NAME)
drugs<-unique(cellLines$DRUG_NAME)
drugsframe<-data.frame(matrix(unlist(drugs), nrow=192, byrow=T))
genesList<-as.list(genes)


cyt<-lapply(genesList, getPValue, geneFile=geneExpression, cellFile=cellLines)


cleanup<- function(name, finalList){
    corrDF<-data.frame(matrix(unlist(finalList), nrow=16224, byrow=T))
    corrDF<-corrDF[,2:3]
    final<-cbind(genes, corrDF)
    names(final)<- c("Gene_Name", "P-Value", "Correlation")
    
    
    final$Correlation<-as.numeric(levels(final$Correlation))[final$Correlation]
    final$`P-Value`<-as.numeric(levels(final$`P-Value`))[final$`P-Value`]
    
    absCorr<-abs(final$Correlation)
    final<-cbind(final,absCorr)
    adjPValue<-p.adjust(final$`P-Value`, method="BH")
    final<-cbind(final, adjPValue)
    fwrite(final, file=name+" Correlation.csv")
}

cleanup(empty2)

Cytarabine<-fread("Cytarabine Correlations.csv")
Vinblastine<-fread("Vinblastine Correlations.csv")
Cisplatin<-fread("Cisplatin Correlations.csv")
temozolomide<-fread("Temozolomide Correlations.csv")
camptothecin<-fread("Camptothecin Correlations.csv")

oncology<-temozolomide %>%
   top_n(-250, adjPValue)%>%
    select(Gene_Name)
fwrite(oncology, "Tmozolomide Oncology.csv")


oncology<-Cytarabine %>%
    top_n(-250, adjPValue)%>%
    select(Gene_Name)
fwrite(oncology, "Cytarabine Oncology.csv")

oncology<-Cisplatin %>%
    top_n(-250, adjPValue)%>%
    select(Gene_Name)
fwrite(oncology, "Cisplatin Oncology.csv")

