library(ggplot2)
library(dplyr)
library(jsonlite)
Camptothecin_Correlationfiltered<-Camptothecin_Correlation%>%
    filter(adjPValue<0.01)%>%
    filter(abs(Correlation)>0.3)
fwrite(Camptothecin_Correlationfiltered, "CamptothecinfilteredCorrelation.csv")

Cytarabine_Correlationfiltered<-Cytarabine_Correlation%>%
    filter(adjPValue<0.01)%>%
    filter(abs(Correlation)>0.3)
fwrite(Cytarabine_Correlationfiltered, "CytarabinefilteredCorrelation.csv")

Temozolomide_Correlationfiltered<-Temozolomide_Correlation%>%
    filter(adjPValue<0.01)%>%
    filter(abs(Correlation)>0.3)
fwrite(Temozolomide_Correlationfiltered, "TemozolomidefilteredCorrelation.csv")

Vinblastine_Correlationfiltered<-Vinblastine_Correlation%>%
    filter(adjPValue<0.01)%>%
    filter(abs(Correlation)>0.3)
fwrite(Vinblastine_Correlationfiltered, "CamptothecinfilteredCorrelation.csv")



Cytarabine_Correlation<-read_csv("Cytarabine Correlation.csv")
Cytarabine_CorrelationCutoff<-Cytarabine_Correlation%>%
    filter(adjPValue<0.05)%>%
    filter(abs(Correlation)>0.3)%>%
    mutate(highorlow=ifelse(Correlation>0, 1, 0))
ggplot(Cytarabine_CorrelationCutoff, aes(adjPValue, Correlation)) +geom_point()
table(Cytarabine_CorrelationCutoff$highorlow)

Temozolomide_Correlation<-read_csv("Temozolomide Correlation.csv")
Temozolomide_CorrelationCutoff<-Temozolomide_Correlation%>%
    filter(adjPValue<0.05)%>%
    filter(abs(Correlation)>0.3)%>%
    mutate(highorlow=ifelse(Correlation>0, 1, 0))
ggplot(Temozolomide_CorrelationCutoff, aes(adjPValue, Correlation)) +geom_point()
table(Temozolomide_CorrelationCutoff$highorlow)

Vinblastine_Correlation<-read_csv("Vinblastine Correlation.csv")
Vinblastine_CorrelationCutoff<-Vinblastine_Correlation%>%
    filter(adjPValue<0.05)%>%
    filter(abs(Correlation)>0.3)%>%
    mutate(highorlow=ifelse(Correlation>0, 1, 0))
ggplot(Vinblastine_CorrelationCutoff, aes(adjPValue, Correlation)) +geom_point()
table(Vinblastine_CorrelationCutoff$highorlow)

VinblastinefilteredCorrelation<-read_csv("VinblastinefilteredCorrelation.csv")
TemozolomidefilteredCorrelation<-read_csv("TemozoleomidefilteredCorrelation.csv")
Cyt<-Cytarabine_Correlationfiltered$Gene_Name
Vin<-Vinblastine_Correlationfiltered$Gene_Name
Camp<-Camptothecin_Correlationfiltered$Gene_Name
Temo<-Temozolomide_Correlationfiltered$Gene_Name
common<-Reduce(intersect, list(Cyt, Vin, Camp, Temo))
common<-as.data.frame(common)
fwrite(common, "Common_genes_0.3.csv")



CytarabinefilteredCorrelationPos<-CytarabinefilteredCorrelation%>%
    filter(Correlation>0)
CamptothecinfilteredCorrelationPos<-CamptothecinfilteredCorrelation%>%
    filter(Correlation>0)
VinblastinefilteredCorrelationPos<-VinblastinefilteredCorrelation%>%
    filter(Correlation>0)
TemozolomidefilteredCorrelationPos<-TemozolomidefilteredCorrelation%>%
    filter(Correlation>0)

CytarabinefilteredCorrelationNeg<-CytarabinefilteredCorrelation%>%
    filter(Correlation<0)
CamptothecinfilteredCorrelationNeg<-CamptothecinfilteredCorrelation%>%
    filter(Correlation<0)
VinblastinefilteredCorrelationNeg<-VinblastinefilteredCorrelation%>%
    filter(Correlation<0)
TemozolomidefilteredCorrelationNeg<-TemozolomidefilteredCorrelation%>%
    filter(Correlation<0)

Cyt<-CytarabinefilteredCorrelationPos$Gene_Name
Vin<-VinblastinefilteredCorrelationPos$Gene_Name
Camp<-CamptothecinfilteredCorrelationPos$Gene_Name
Temo<-TemozolomidefilteredCorrelationPos$Gene_Name
common<-Reduce(intersect, list(Cyt, Vin, Camp, Temo))
common<-as.data.frame(common)
fwrite(common, "Common Pos Correlation genes.csv")

Cyt<-CytarabinefilteredCorrelationNeg$Gene_Name
Vin<-VinblastinefilteredCorrelationNeg$Gene_Name
Camp<-CamptothecinfilteredCorrelationNeg$Gene_Name
Temo<-TemozolomidefilteredCorrelationNeg$Gene_Name
common<-Reduce(intersect, list(Cyt, Vin, Camp, Temo))
common<-as.data.frame(common)
fwrite(common, "Common Neg Correlation genes.csv")



GOAnalysis<-fromJSON("GO_Common_0.3.json")%>%as.data.frame
GOAnalysisFinal<-cbind(GOAnalysis2$overrepresentation.result.term$label,
                       GOAnalysis2$overrepresentation.result.input_list$number_in_list,
                       GOAnalysis2$overrepresentation.result.input_list$fold_enrichment, 
                       GOAnalysis2$overrepresentation.result.input_list$fdr, 
                       GOAnalysis2$overrepresentation.result.input_list$expected,
                       GOAnalysis2$overrepresentation.result.input_list$pValue)
colnames(GOAnalysisFinal)<-c( "Function","Actual", "Fold Enrichment", "FDR", "Expected", "PValue")
fwrite(GOAnalysisFinal, "GOAnalysis_0.3.csv")
