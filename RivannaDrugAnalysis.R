install.packages("readr")
install.packages

library(dplyr)
library(readr)
library(wPerm)
library(data.table)
library(parallel)
library(doParallel)


getAllDrugs<-function(){
    cellLines<-readxl::read_xlsx("GDSC2_fitted_dose_response_25Feb20.xlsx") #Read cell line file
    cellLines<-cellLines%>%   #Get only relevant columns
        select("CELL_LINE_NAME", "DRUG_NAME", "LN_IC50")
    geneExpression<-read_tsv("CosmicCLP_CompleteGeneExpression.tsv")  #Read gene expression file
    geneExpression<-geneExpression%>%  #Get only relevant columns
        select(-"SAMPLE_ID", -"REGULATION")
    genes<-unique(geneExpression$GENE_NAME)  #Get a list of unique genes
    drugs<-unique(cellLines$DRUG_NAME)  #Get a list of unique drugs
    genesList<-as.list(genes)  #Ensure genes is a list
    drugsList<-as.list(drugs)
    mclapply(drugsList, drugCorTest, geneFile=geneExpression, cellFile=cellLines, mc.cores = 7)  #Apply DrugcorTest to all drugs
}

drugCorTest<-function(geneFile, cellFile, drug){ #For single drug
    cellFile<-cellLines
    name<-"Camptothecin"
    name<-drug  #Get name of drug (for future use and file naming)
    selectLines<-cellFile %>%  #Get cell lines treated with that drug
        filter(DRUG_NAME == name)
    raw2<-lapply(genesList, getCorr, geneFile=geneExpression, selectLines=selectLines) #Apply getCorr to all genes
    cleanup(name, raw)
}

getCorr<-function(geneFile, selectLines, geneName){ #Get correlation between the given cell lines and all genes
    #geneName<-genesList[1]
    #geneFile<-geneExpression
    
    selectGenes<-geneFile %>%  #Get all genes of one name
        filter(GENE_NAME== geneName)
    merged<-left_join(selectLines,selectGenes, by=c("CELL_LINE_NAME"="SAMPLE_NAME"))  #Merge the cell line and gene files
    merged<-na.omit(merged)  #Remove all NA values (Cell lines without the given gene)
    x<-merged$LN_IC50 #Set x as the IC50 values
    y<-merged$Z_SCORE  #Set y as the z score expression values
    cor.test(x,y, method="spearman")  #Do correlation test
}
cleanup<- function(name, finalList){
    finalList<-raw2
    corrDF<-data.frame(matrix(unlist(finalList), nrow=length(finalList), byrow=T))  
    corrDF<-corrDF[,2:3]
    final<-cbind(genes, corrDF)
    names(final)<- c("Gene_Name", "P-Value", "Correlation")
    final$Correlation<-as.numeric(final$Correlation)
    final$`P-Value`<-as.numeric(final$`P-Value`)
    #final$Correlation<-as.numeric(levels(final$Correlation))[final$Correlation]
    #final$`P-Value`<-as.numeric(levels(final$`P-Value`))[final$`P-Value`]
    
    absCorr<-abs(final$Correlation)
    final<-cbind(final,absCorr)
    adjPValue<-p.adjust(final$`P-Value`, method="BH")
    final<-cbind(final, adjPValue)
    fwrite(final, file="Camptothecin Correlation.csv")
}

drugCorTest(geneExpression, cellLines, "Camptothecin")
