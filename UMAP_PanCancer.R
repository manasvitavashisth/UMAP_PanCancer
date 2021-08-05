rm(list=ls())
library(Seurat)
library(dplyr)
library(stringr)
library(ggrepel)
library(survminer)
library(survival)

cancerListNames = list(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  #"Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)
cancerListNames1 = c(
  "Liver",
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  #"Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

cancerListNames=list(
  "Thymoma",
  "Lung_Adenocarcinoma",
  "Breast",
  "Lower_Grade_Glioma",
  "Sarcoma",
  "Pancreatic",
  "Prostate",
  "Kidney_Chromophobe",
  "Adrenocortical",
  "Bile_Duct",
  "Glioblastoma",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Thyroid",
  "Endometrioid",
  "Bladder",
  "Ocular_Melanomas",
  "Large_B-cell_Lymphoma",
  "Esophageal",
  "Head_and_Neck",
  "Melanoma",
  "Mesothelioma",
  "Ovarian",
  "Rectal",
  "Cervical",
  "Colon",
  "Lung_Squamous_Cell",
  "Pheochromocytoma_Paraganglioma",
  "Stomach",
  "Testicular",
  "Uterine"
)

cancerListNames1=c(
  "Liver",
  "Thymoma",
  "Lung_Adenocarcinoma",
  "Breast",
  "Lower_Grade_Glioma",
  "Sarcoma",
  "Pancreatic",
  "Prostate",
  "Kidney_Chromophobe",
  "Adrenocortical",
  "Bile_Duct",
  "Glioblastoma",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Thyroid",
  "Endometrioid",
  "Bladder",
  "Ocular_Melanomas",
  "Large_B-cell_Lymphoma",
  "Esophageal",
  "Head_and_Neck",
  "Melanoma",
  "Mesothelioma",
  "Ovarian",
  "Rectal",
  "Cervical",
  "Colon",
  "Lung_Squamous_Cell",
  "Pheochromocytoma_Paraganglioma",
  "Stomach",
  "Testicular",
  "Uterine"
)

cancerListNames=list(
  "Thymoma",
  "Lung_Adenocarcinoma",
  "Breast",
  "Lower_Grade_Glioma",
  "Sarcoma",
  "Pancreatic",
  "Prostate",
  "Kidney_Chromophobe",
  "Adrenocortical",
  "Bile_Duct",
  "Glioblastoma",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Thyroid",
  "Endometrioid",
  "Bladder",
  "Ocular_Melanomas",
  "Large_B-cell_Lymphoma",
  "Esophageal",
  "Head_and_Neck",
  "Melanoma",
  "Mesothelioma",
  "Ovarian",
  "Rectal",
  "Cervical",
  "Colon",
  "Lung_Squamous_Cell",
  "Pheochromocytoma_Paraganglioma",
  "Stomach",
  "Testicular",
  "Uterine"
)

cancerListNames1=c(
  "Liver",
  "Thymoma",
  "Lung_Adenocarcinoma",
  "Breast",
  "Lower_Grade_Glioma",
  "Sarcoma",
  "Pancreatic",
  "Prostate",
  "Kidney_Chromophobe",
  "Adrenocortical",
  "Bile_Duct",
  "Glioblastoma",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Thyroid",
  "Endometrioid",
  "Bladder",
  "Ocular_Melanomas",
  "Large_B-cell_Lymphoma",
  "Esophageal",
  "Head_and_Neck",
  "Melanoma",
  "Mesothelioma",
  "Ovarian",
  "Rectal",
  "Cervical",
  "Colon",
  "Lung_Squamous_Cell",
  "Pheochromocytoma_Paraganglioma",
  "Stomach",
  "Testicular",
  "Uterine"
)

cancerListNames=list(
  "Liver",
  "Thymoma",
  "Lung_Adenocarcinoma",
  "Breast",
  "Lower_Grade_Glioma",
  "Sarcoma",
  "Pancreatic",
  "Prostate",
  "Kidney_Chromophobe",
  "Adrenocortical",
  "Bile_Duct",
  "Glioblastoma",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Thyroid",
  "Endometrioid",
  "Bladder",
  
  "Large_B-cell_Lymphoma",
  
  "Testicular"
)

setwd("/Applications/UPenn/Summer2020/TCGA_Data")
plotdf=read.table("MedianSurvival", header = TRUE,sep=" ",stringsAsFactors = FALSE,na.strings=c("", "NA"))
maxoverlap=read.table("Tumor_GeneOverlap", header = TRUE,sep=" ",stringsAsFactors = FALSE,na.strings=c("", "NA"))
maxoverlap[,2]=str_replace_all(maxoverlap[,2],'[+]','')
maxoverlap[,2]=str_replace_all(maxoverlap[,2],'[ ]','')
plotdf$Cancer=str_replace_all(plotdf$Cancer,'[_]','')
FOXM1fit1$cancername=str_replace_all(FOXM1fit1$V3,'[_]','')

plotdf1=plotdf[plotdf$Cancer!="PheochromocytomaParaganglioma"& plotdf$Cancer!="Testicular" &plotdf$Cancer!="LargeB-cellLymphoma" & plotdf$Cancer!="Thymoma",]


avgslope=mean(FOXM1fit[FOXM1fit$V2>0.5,"V1"])
stdevavgslope=sd(FOXM1fit[FOXM1fit$V2>0.5,"V1"])
FOXM1fit$deviation=FOXM1fit$V1>(avgslope+stdevavgslope) | FOXM1fit$V1<(avgslope-stdevavgslope)

plotdf$exclusion=plotdf$Cancer %in% exclusionlist_rsq
plotdf$goodfit=plotdf$Rsq >0.5
t.test(plotdf[plotdf$goodfit,"median_survival1"],plotdf[plotdf$goodfit==FALSE,"median_survival1"])

boxplot(median_survival~goodfit,data=plotdf, main="Median OS",
        xlab="Rsq for Forced Fit >0.5", ylab="Avg Median OS of tumor cohort(days)")


plotdf1$exclusion=plotdf1$Cancer %in% exclusionlist_rsq
plotdf1$goodfit=plotdf1$Rsq >0.5

plotdf2=plotdf1[plotdf$median_survival_zeros<10000,]

t.test(plotdf2[plotdf2$exclusion,"median_survival"],plotdf2[plotdf2$exclusion==FALSE,"median_survival"])


t.test(plotdf1[plotdf1$goodfit,"median_survival1"],plotdf1[plotdf1$goodfit==FALSE,"median_survival1"])

boxplot(geneoverlap_survival1[geneoverlap_survival1$exclusion_rsq,"median_survival"],geneoverlap_survival1[geneoverlap_survival1$exclusion_rsq==FALSE,"median_survival"])

boxplot(median_survival1~goodfit,data=plotdf1, main="Median OS",
        xlab="Rsq>0.5", ylab="Avg Median OS of tumor cohort(years)")


fitslope=merge(plotdf,FOXM1fit1,by.x = "Cancer",by.y = "V3")
fitslope1=fitslope[fitslope$median_survival_zeros>-6000,]

ggplot(fitslope,aes(V2,median_survival,label=Cancer))+geom_point(color = "red")+geom_label_repel()+xlab("Fitted Slope Rsq")
ggplot(fitslope,aes(V2,median_survival,label=Cancer))+geom_point(color = "red")+xlab("Fitted Slope Rsq")+geom_smooth(method = "lm", se = FALSE)

ggplot(fitslope1,aes(V2,median_survival,label=Cancer))+geom_point(color = "red")+geom_label_repel()+xlab("Fitted Slope Rsq")
ggplot(fitslope1,aes(V2,median_survival,label=Cancer))+geom_point(color = "red")+xlab("Fitted Slope Rsq")+geom_smooth(method = "lm", se = FALSE)

fitslope1$goodfit=fitslope1$V2>0.5
boxplot(median_survival~badfit,data=fitslope1, main="Median Survival",
        xlab="Good fit for fitted slope (Rsq>0.5)", ylab="Avg Median survival of tumor cohort(days)")



exclusionlist_overlap=str_replace_all(exclusionlist_overlap,'[_]','')
exclusionlist_rsq=str_replace_all(exclusionlist_rsq,'[_]','')
FOXM1fit$tumor=str_replace_all(FOXM1fit$V3,'[_]','')

geneoverlap_survival=merge(maxoverlap,plotdf,by.x="Added.Cancer",by.y="Cancer")

rsq=merge(geneoverlap_survival1,FOXM1fit,by.x="Added.Cancer",by.y="tumor")

ggplot(rsq,aes(V1,median_survival_zeros,label=Added.Cancer))+geom_point(color = "red")+geom_label_repel()

ggplot(geneoverlap_survival1,aes(log(Number.of.Genes.In.Common),median_survival_zeros,label=Added.Cancer))+geom_point(color = "red")+geom_label_repel()

ggplot(geneoverlap_survival,aes(log(Number.of.Genes.In.Common),median_survival_zeros,label=Added.Cancer))+geom_point(color = "red")+geom_smooth(method = "lm", se = FALSE)

ggplot(geneoverlap_survival1,aes(log(Number.of.Genes.In.Common),median_survival_zeros,label=Added.Cancer))+geom_point(color = "red")+geom_smooth(method = "lm", se = FALSE)

ggplot(geneoverlap_survival1,aes(Rsq_FittedSlope,median_survival_zeros,label=Added.Cancer))+geom_point(color = "red")+geom_label_repel()

ggplot(geneoverlap_survival1,aes(Rsq_FittedSlope,median_survival_zeros,label=Added.Cancer))+geom_point(color = "red")+geom_smooth(method = "lm", se = FALSE)


geneoverlap_survival1=geneoverlap_survival[geneoverlap_survival$median_survival_zeros>-6000,]
geneoverlap_survival1=geneoverlap_survival1[geneoverlap_survival1$Added.Cancer != "PheochromocytomaParaganglioma",]
geneoverlap_survival$exclusion=geneoverlap_survival$Added.Cancer %in% exclusionlist_overlap

geneoverlap_survival1$exclusion=geneoverlap_survival1$Rsq_FittedSlope<0.5
geneoverlap_survival1$exclusion_rsq=geneoverlap_survival1$Rsq<0.5


geneoverlap_survival1$exclusion_rsq=geneoverlap_survival1$Added.Cancer %in% exclusionlist_overlap

avgslope=mean(geneoverlap_survival[geneoverlap_survival$exclusion==FALSE & geneoverlap_survival$Rsq>0.5,"Exponent"])
stdevavgslope=sd(geneoverlap_survival[geneoverlap_survival$exclusion==FALSE & geneoverlap_survival$Rsq>0.5,"Exponent"])

mean(geneoverlap_survival1[geneoverlap_survival1$exclusion==FALSE,"median_survival"])

t.test(geneoverlap_survival1[geneoverlap_survival1$exclusion_rsq,"median_survival"],geneoverlap_survival1[geneoverlap_survival1$exclusion_rsq==FALSE,"median_survival"])

t.test(rsq[rsq$deviation,"median_survival"],rsq[rsq$deviation==FALSE,"median_survival"])

boxplot(geneoverlap_survival1[geneoverlap_survival1$exclusion_rsq,"median_survival"],geneoverlap_survival1[geneoverlap_survival1$exclusion_rsq==FALSE,"median_survival"])

boxplot(median_survival~exclusion_rsq,data=geneoverlap_survival1, main="Median PFI",
        xlab="Overlap exclusion", ylab="Avg Median PFI of tumor cohort(days)")

boxplot(median_survival~deviation,data=rsq, main="Median Survival",
        xlab="Deviation from Mean Exponent", ylab="Avg Median survival of tumor cohort(days)")

CNN1=read.table("trak1_Liver_LMNB1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))
# df.map = read.table("Data_Liver_HiSeqV2", header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
# df.map1=read.table("Data_Lung_Adenocarcinoma_HiSeqV2", header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
# df.map=cbind(df.map,df.map1[,-1])
CNN1=CNN1[CNN1[,5] %in% all.genes,]

setwd("/Applications/UPenn/Summer2020/TCGA_Data/TCGA_data")
df.map1 = read.table("Data_Liver_HiSeqV2", header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
genes=df.map1[,1]
df.map=df.map1[,substr(colnames(df.map1),14,15)=='11']
patient=dim(df.map)[2]
a=substr(colnames(df.map),1,12)
df.map2=df.map1[,(substr(colnames(df.map1),1,12) %in% a) & (substr(colnames(df.map1),14,15)=='01')]
df.map=cbind(df.map,df.map2)
patient[2]=dim(df.map2)[2]
df.surv = read.table("Data_Liver_survival.txt", header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
trak1CancerListNames1 = paste("Data_", cancerListNames, "_survival.txt", sep = "")
df.surv=df.surv[df.surv$sample %in% a,]
#df.residual=(read.table("ResidualswrtFittedSlope_Liver_FOXM1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA")))
#df.residual1=read.table("Residuals_Liver_FOXM1.txt", header = TRUE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))

fit = survfit(Surv(df.surv$OS.time, df.surv$OS) ~ 1, data = df.surv)
ggsurvplot(
  survfit(Surv(df.surv$PFI.time, df.surv$PFI) ~ 1, data = df.surv),
  xlab = "Days",
  ylab = "Progression free interval probability")
median_survival =surv_median(fit)$median
median_survival_zeros =surv_median(fit)$median
#cancerListNames = list("Lung_Adenocarcinoma","Breast","Bladder","Thyroid","Uterine","Stomach")
OS <- data.frame(matrix(ncol = 5, nrow = 1))
PFI <- data.frame(matrix(ncol = 5, nrow = 1))
k=1
m=1
for(i in 1:dim(df.map)[2])
{
  for(j in 1:dim(df.surv)[1])
  {
    if(a[i]==df.surv$sample[j] & df.surv$OS[j]!=0)
    {
      OS[k,1]=df.surv$OS.time[j]
      OS[k,2]=df.residual[i]
      OS[k,3]=df.surv$sample[j]
      OS[k,4]="Liver"
      OS[k,5]=df.residual1$Studentized[i]
      k=k+1
      break
    }
    if(a[i]==df.surv$sample[j] & df.surv$PFI[j]!=0)
    {
      PFI[m,1]=df.surv$PFI.time[j]
      PFI[m,2]=df.residual[i]
      PFI[m,3]=df.surv$sample[j]
      PFI[m,4]="Liver"
      PFI[m,5]=df.residual1$Studentized[i]
      m=m+1
      break
    }
  }
}

plot(PFI[,1],abs(PFI[,2]))
fit=lm(X1~abs(X5),data=PFI)
# Trk[,2]=normalize(Trk[,1], method = "standardize", range = c(0, 1), margin = 2)
# surv_all=Trk
trak1CancerListNames = paste("Data_", cancerListNames, "_HiSeqV2", sep = "")
df.residual=(read.table("ResidualswrtFittedSlope_Liver_FOXM1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA")))

residualCancerListNames = paste("ResidualswrtFittedSlope_", cancerListNames, "_FOXM1.txt", sep = "")
#trak1CancerListNames1 = paste("Data_", cancerListNames, "_survival.txt", sep = "")

k=3
for (i in 1:length(cancerListNames)) {
  df.map1 = read.table(trak1CancerListNames[i], header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
  df.map3=df.map1[,substr(colnames(df.map1),14,15)=='11']
  #a=str_replace_all(colnames(df.map1),'[.]','-')
  patient[k]=dim(df.map3)[2]
  k=k+1
  df.map=cbind(df.map,df.map3)
  #patient=dim(df.map)[2]
  a=substr(colnames(df.map3),1,12)
  df.map2=df.map1[,(substr(colnames(df.map1),1,12) %in% a) & (substr(colnames(df.map1),14,15)=='01')]
  df.map=cbind(df.map,df.map2)
  patient[k]=dim(df.map2)[2]
  k=k+1
  # df.residual1=read.table(residualCancerListNames[i], header = FALSE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))
  # df.residual=cbind(df.residual,df.residual1)
  # df.surv1=read.table(trak1CancerListNames1[i], header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
  # df.surv1=df.surv1[df.surv1$sample %in% a,]
  # fit = survfit(Surv(df.surv1$OS.time, df.surv1$OS) ~ 1, data = df.surv1)
  # median_survival[i+1] =surv_median(fit)$median
  # median_survival_zeros[i+1] =surv_median(fit)$median
  # if(is.na(surv_median(fit)$median))
  # {
  #   median_survival[i+1]=max(fit$time)*fit$surv[length(fit$surv)]/0.5
  #   median_survival_zeros[i+1] =100000
  # }
  # df.surv = read.table(trak1CancerListNames1[i], header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
  # Trk <- data.frame(matrix(ncol = 2, nrow = 1))
  # for(i in 1:dim(df.map1)[2])
  # {
  #   for(j in 1:dim(df.surv)[1])
  #   {
  #     if(substr(colnames(df.map1)[i],6,7)==substr(df.surv$X_PATIENT[j],6,7) & substr(colnames(df.map1)[i],9,12)==substr(df.surv$X_PATIENT[j],9,12))
  #     {
  #       Trk[i,1]=df.surv$OS.time[j]
  #       break
  #     }
  #   }
  # }
  # Trk[,2]=normalize(Trk[,1], method = "standardize", range = c(0, 1), margin = 2)
  # surv_all=rbind(surv_all,Trk)
}



setwd("/Applications/UPenn/Summer2020/TCGA_Data/TCGA_data")
df.map1 = read.table("Data_Liver_HiSeqV2", header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
genes=df.map1[,1]
df.map=df.map1[,substr(colnames(df.map1),14,15)=='01']
patient=dim(df.map)[2]
trak1CancerListNames = paste("Data_", cancerListNames, "_HiSeqV2", sep = "")
for (i in 1:length(cancerListNames)) {
  df.map1 = read.table(trak1CancerListNames[i], header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
  df.map3=df.map1[,substr(colnames(df.map1),14,15)=='01']
  #a=str_replace_all(colnames(df.map1),'[.]','-')
  patient[i+1]=dim(df.map3)[2]
  df.map=cbind(df.map,df.map3)
  
}

trak1CancerListNames = paste("Data_", cancerListNames4, "_HiSeqV2", sep = "")
for (i in 1:length(cancerListNames4)) {
  df.map1 = read.table(trak1CancerListNames[i], header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
  df.map3=df.map1[,substr(colnames(df.map1),14,15)=='11']
  #a=str_replace_all(colnames(df.map1),'[.]','-')
  patient[i+32]=dim(df.map3)[2]
  df.map=cbind(df.map,df.map3)
  
}

FOXM1fit1=read.table("Intercept_RsqwrtFittedSlope.txt", header = FALSE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))

FOXM1fit=read.table("FOXM1fit_slope_rsq.txt", header = FALSE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))
colnames(FOXM1fit)=c("Exponent","Rsq","Cancer")
colnames(FOXM1fit1)=c("Intercept","Rsq_FittedSlope","CancerNames")

median_survival1=median_survival/365
median_survival_zeros1=median_survival_zeros/365
cancerListPatients=rep((cancerListNames1),times=as.vector(patient))
FOXM1fitrep=rep((FOXM1fit1[,2]),times=as.vector(patient))
FOXM1fitrep_exp=rep((FOXM1fit[,2]),times=as.vector(patient))
surv_median_label=rep((median_survival),times=as.vector(patient))
surv_median_zeros_label=rep((median_survival_zeros),times=as.vector(patient))
surv_median_label=rep((median_survival1),times=as.vector(patient))
surv_median_zeros_label=rep((median_survival_zeros1),times=as.vector(patient))

 
trak1CancerListNames1 = paste("kmsig2_", cancerListNames, "_coxph", sep = "")
df.coxph=read.table("kmsig2_Liver_coxph", header = TRUE,sep=" ",stringsAsFactors = FALSE,na.strings=c("", "NA"))
FOXM1=cbind(df.coxph[df.coxph$Gene=="FOXM1",],cancerListNames1[1])
TOP2A=cbind(df.coxph[df.coxph$Gene=="TOP2A",],cancerListNames1[1])
LMNB1=cbind(df.coxph[df.coxph$Gene=="LMNB1",],cancerListNames1[1])
for (i in 1:length(cancerListNames)){
  df.coxph=read.table(trak1CancerListNames1[i], header = TRUE,sep=" ",stringsAsFactors = FALSE,na.strings=c("", "NA"))
  FOXM1[i+1,]=cbind(df.coxph[df.coxph$Gene=="FOXM1",],cancerListNames1[i+1])
  TOP2A[i+1,]=cbind(df.coxph[df.coxph$Gene=="TOP2A",],cancerListNames1[i+1])
  LMNB1[i+1,]=cbind(df.coxph[df.coxph$Gene=="LMNB1",],cancerListNames1[i+1])
}
FOXM1_sig=FOXM1[FOXM1$p.value_Coeff<0.05,]
TOP2A_sig=TOP2A[TOP2A$p.value_Coeff<0.05,]
LMNB1_sig=LMNB1[LMNB1$p.value_Coeff<0.05,]

foxm1_overlap=FOXM1[FOXM1$`cancerListNames1[1]` %notin% exclusionlist_overlap,]
top2a_overlap=TOP2A[TOP2A$`cancerListNames1[1]` %notin% exclusionlist_overlap,]
LMNB1_overlap=LMNB1[LMNB1$`cancerListNames1[1]` %notin% exclusionlist_overlap,]
a=LMNB1_overlap[LMNB1_overlap$p.value_Coeff<0.05,]
b=top2a_overlap[top2a_overlap$p.value_Coeff<0.05,]
c=foxm1_overlap[foxm1_overlap$p.value_Coeff<0.05,]

library('plotrix')

pie(c(7,1,9),radius=0.2,init.angle=90)

floating.pie(0,0,c(9,8,0),radius=0.8,startpos =pi/2,col=c("firebrick1","slategray1"))
floating.pie(0,0,c(8,8,1),radius=0.6,startpos =pi/2,col=c("firebrick1","slategray1","lightseagreen"))
floating.pie(0,0,c(7,9,1),radius=0.4,startpos =pi/2,col=c("firebrick1","slategray1","lightseagreen"))

floating.pie(0,0,c(9,8,0),radius=0.8,startpos =pi/2,col=c("darkturquoise","slategray1"))
floating.pie(0,0,c(8,8,1),radius=0.6,startpos =pi/2,col=c("darkturquoise","slategray1","midnightblue"))
floating.pie(0,0,c(7,9,1),radius=0.4,startpos =pi/2,col=c("darkturquoise","slategray1","midnightblue"))


floating.pie(0,0,c(9,0,8),radius=0.8,startpos =pi/2)
floating.pie(0,0,c(8,1,8),radius=0.6,startpos =pi/2)
floating.pie(0,0,c(7,1,9),radius=0.4,startpos =pi/2)

floating.pie(0,0,c(15,1),radius=0.8,startpos =pi/2,col=c("slategray1","midnightblue"))
floating.pie(0,0,c(15,1),radius=0.6,startpos =pi/2,col=c("slategray1","midnightblue"))
floating.pie(0,0,c(15,1),radius=0.4,startpos =pi/2,col=c("slategray1","midnightblue"))

floating.pie(0,0,c(0,14),radius=1.4,startpos =pi/2,col=c("aquamarine2","white"))
floating.pie(0,0,c(1,14),radius=1.5,startpos =pi/2,col=c("aquamarine2","white"))

floating.pie(0,0,c(2,13),radius=0.3*pi,startpos =0,col=c("aquamarine2","white"))
floating.pie(0,0,c(1,14),radius=0.2*pi,startpos =2*pi/15,col=c("aquamarine2","white"))
floating.pie(0,0,c(1,14),radius=0.1*pi,startpos =0,col=c("aquamarine2","white"))

floating.pie(0,0,c(2,13),radius=0.2*pi,startpos =pi/2,col=c("aquamarine2","white"))
floating.pie(0,0,c(1,14),radius=0.3*pi,startpos =pi/2,col=c("aquamarine2","white"))
floating.pie(0,0,c(1,14),radius=0.4*pi,startpos =pi/2,col=c("aquamarine2","white"))



`%notin%`=Negate(`%in%`)
FOXM1=subset(df.coxph,df.coxph$Gene=="FOXM1")
df.residual=read.table("Residuals_Liver_KIF18B.txt", header = TRUE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))
residualCancerListNames = paste("Residuals_", cancerListNames, "_KIF18B.txt", sep = "")
for (i in 1:length(cancerListNames)){
df.residual1=read.table(residualCancerListNames[i], header = TRUE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))
df.residual=rbind(df.residual,df.residual1)
}
# write.table(surv_all, file = "SurvTime_InOrderofAllCancers_Normalized", append = FALSE, quote = TRUE, sep = " ",
#             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#             col.names = TRUE, qmethod = c("escape", "double"),
#             fileEncoding = "")
# 

df.surv = read.table("Data_Liver_survival.txt", header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
trak1CancerListNames1 = paste("Data_", cancerListNames, "_survival.txt", sep = "")
fit = survfit(Surv(df.surv$OS.time, df.surv$OS) ~ 1, data = df.surv)
ggsurvplot(
  survfit(Surv(df.surv$OS.time, df.surv$OS) ~ 1, data = df.surv),
  xlab = "Days",
  ylab = "Overall survival probability")
median_survival =surv_median(fit)$median
median_survival_zeros =surv_median(fit)$median
for (i in 1:length(cancerListNames)){
  df.surv1=read.table(trak1CancerListNames1[i], header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
  fit = survfit(Surv(df.surv1$OS.time, df.surv1$OS) ~ 1, data = df.surv1)
  median_survival[i+1] =surv_median(fit)$median
  median_survival_zeros[i+1] =surv_median(fit)$median
  if(is.na(surv_median(fit)$median))
  {
    median_survival[i+1]=max(fit$time)
    median_survival_zeros[i+1] =-6000
  }
}

FOXM1fit1=read.table("Intercept_RsqwrtFittedSlope.txt", header = FALSE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))

FOXM1fit=read.table("FOXM1fit_slope_rsq.txt", header = FALSE,sep=",",stringsAsFactors = FALSE,na.strings=c("", "NA"))

median_survival1=median_survival/365
median_survival_zeros1=median_survival_zeros/365

avgslope=mean(FOXM1fit[FOXM1fit$V2>0.5,"V1"])
stdevavgslope=sd(FOXM1fit[FOXM1fit$V2>0.5,"V1"])
FOXM1fit$deviation=FOXM1fit$V1>(avgslope+stdevavgslope) | FOXM1fit$V1<(avgslope-stdevavgslope)
  

colnames(FOXM1fit)=c("Exponent","Rsq","Cancer")
colnames(FOXM1fit1)=c("Intercept","Rsq_FittedSlope","CancerNames")

plot(FOXM1fit$Rsq,median_survival)

plot(FOXM1fit$Rsq,median_survival_zeros)

plot(FOXM1fit1$Rsq_FittedSlope,median_survival)

plot(FOXM1fit1$Rsq_FittedSlope,median_survival_zeros)
text(FOXM1fit1$Rsq_FittedSlope,median_survival_zeros,FOXM1fit1$Cancer)
ggplot(df.input,aes(V3,V1,label=V5))+geom_point(color = "blue")+geom_label_repel()
plotdf=cbind(FOXM1fit,FOXM1fit1,median_survival1,median_survival_zeros1)
ggplot(plotdf,aes(Rsq_FittedSlope,median_survival,label=CancerNames))+geom_point(color = "red")+geom_label_repel()
ggplot(plotdf,aes(Rsq_FittedSlope,median_survival,label=CancerNames))+geom_point(color = "red")+geom_smooth(method = "lm", se = FALSE)

fit=lm(plotdf$median_survival~plotdf$Rsq_FittedSlope,data=plotdf)

plotdf1=plotdf[plotdf$Cancer!="PheochromocytomaParaganglioma"& plotdf$Cancer!="Testicular" &plotdf$Cancer!="LargeB-cellLymphoma" & plotdf$Cancer!="Thymoma",]

plotdf1=plotdf[plotdf$median_survival_zeros>-6000,]
ggplot(plotdf1,aes(Rsq_FittedSlope,median_survival,label=CancerNames))+geom_point(color = "red")+geom_label_repel()
ggplot(plotdf1,aes(Rsq_FittedSlope,median_survival,label=CancerNames))+geom_point(color = "red")+geom_smooth(method = "lm", se = FALSE)

ggplot(plotdf1,aes(Rsq,median_survival_zeros,label=CancerNames))+geom_point(color = "red")+geom_smooth(method = "lm", se = FALSE)
fit=lm(plotdf1$median_survival_zeros~plotdf1$Rsq,data=plotdf1)
fitCancerListNames1 = paste("_", cancerListNames1, "_", sep = "")
FOXM1fit=merge(fitCancerListNames1,FOXM1fit1,by.x=fitCancerListNames1,by.y=FOXM1fit1$V3)

write.table(plotdf, file = "PFIMedianSurvival", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

for (j in 1:31)
{
  for (i in 1:31)
  {
    if(trak1CancerListNames1[j]==fitCancerListNames1[i+1])
    {
      FOXM1fit[i+1,]=FOXM1fit1[j,]
    }
  }
}
cancerListPatients=rep((cancerListNames1),times=as.vector(patient))

cancerListPatients=rep((cancerListNames2),times=as.vector(patient))
x=rep(c("Adjacent","Tumor"),length(cancerListNames1))
cancerListPatients2=rep(x,times=as.vector(patient))

cancerListPatients2=rep(c("Tumor","Adjacent"),times=c(9112,655))

cancerListNames2 = c(
  #"Acute_Myeloid_Leukemia" #Problematic, no primary tumor site
  #"Adrenocortical",
  "Liver",
  "Tumor Liver",
  "Bile Duct",
  "Tumor Bile Duct",
  "Bladder",
  "Tumor Bladder",
  "Breast",
  "Tumor Breast",
  #"Cervical",
  #"Colon_and_Rectal",
  "Colon",
  "Tumor Colon",
  "Endometrioid",
  "Tumor Endometrioid",
  "Esophageal",
  "Tumor Esophageal",
  #"Glioblastoma",
  "Head and Neck",
  "Tumor Head and Neck",
  "Kidney Chromophobe",
  "Tumor Kidney Chromophobe",
  "Kidney Clear Cell",
  "Tumor Kidney Clear Cell",
  "Kidney Papillary Cell",
  "Tumor Kidney Papillary Cell",
  #"Large_B-cell_Lymphoma",
  #"Liver",
  #"Lower_Grade_Glioma_and_Glioblastoma",
  #"Lower_Grade_Glioma",
  "Lung Adenocarcinoma",
  "Tumor Lung Adenocarcinoma",
  #"Lung",
  "Lung Squamous Cell",
  "Tumor Lung Squamous Cell",
  #"Melanoma",
  #"Mesothelioma",
  #"Ocular_Melanomas",
  #"Ovarian",
  #"Pancreatic",
  #"Pheochromocytoma_Paraganglioma",
  #"Prostate",
  "Rectal",
  "Tumor Rectal",
  #"Sarcoma",
  "Stomach",
  "Tumor Stomach",
  #"Testicular",
  #"Thymoma",
  "Thyroid",
  "Tumor Thyroid"
  #"Uterine"
)


FOXM1fitrep=rep((FOXM1fit1[,2]),times=as.vector(patient))
FOXM1fitrep_exp=rep((FOXM1fit[,2]),times=as.vector(patient))
surv_median_label=rep((median_survival),times=as.vector(patient))
surv_median_zeros_label=rep((median_survival_zeros),times=as.vector(patient))
surv_median_label=rep((median_survival1),times=as.vector(patient))
surv_median_zeros_label=rep((median_survival_zeros1),times=as.vector(patient))

surv_median_label1=rep((plotdf$median_survival),times=as.vector(patient))
surv_median_zeros_label1=rep((plotdf$median_survival_zeros),times=as.vector(patient))


all.genes <- rownames(pbmc)
all.genes=CNN1$V5
pbmc <- RunUMAP(pbmc,features = rownames(pbmc))
pbmc <- RunUMAP(pbmc,features = CNN1)

DimPlot(pbmc, reduction = "umap",group.by = "cancer",label=TRUE,repel = TRUE,label.size = 5)

exclusionlist_Liver=c("Large_B-cell_Lymphoma","Head_and_Neck","Testicular","Stomach","Mesothelioma","Rectal","Melanoma","Esophageal","Ovarian","Cervical","Colon","Lung_Squamous_Cell","Pheochromocytoma_Paraganglioma","Uterine")
exclusionlabel=cancerListPatients %in% exclusionlist_Liver

exclusionlist_overlap=c("Ocular_Melanomas","Large_B-cell_Lymphoma","Head_and_Neck","Testicular","Stomach","Mesothelioma","Rectal","Melanoma","Esophageal","Ovarian","Cervical","Colon","Lung_Squamous_Cell","Pheochromocytoma_Paraganglioma","Uterine")

exclusionlabel=cancerListPatients %in% exclusionlist_overlap

exclusionlist_rsq=c("Stomach","Colon","Rectal","Esophageal","PheochromocytomaParaganglioma","OcularMelanomas","LungSquamous_Cell","HeadandNeck","Mesothelioma","Ovarian","Uterine","Melanoma","Cervical")

exclusionlist_rsq=c("Stomach","Colon","Rectal","Esophageal","Pheochromocytoma_Paraganglioma","Ocular_Melanomas","Lung_Squamous_Cell","Head_and_Neck","Mesothelioma","Ovarian","Uterine","Melanoma","Cervical")


exclusionlabelrsq=cancerListPatients %in% exclusionlist_rsq
exclusionlabelrsq1=cancerListPatients %in% exclusionlist_rsq

exclusionlabelrsq1[exclusionlabelrsq]="Dysregulated"
exclusionlabelrsq1[!exclusionlabelrsq]="Regulated"

#median_survival_zeros[median_survival_zeros==1000000]=10000

CNN1=as.vector(CNN1)
CNN1=t(CNN1)
#genes=df.map[,1]
#df.map=df.map[,substr(colnames(df.map),14,15)=='01']
rownames(df.map)=genes
pbmc <- CreateSeuratObject(counts = df.map)
pbmc[["cancer"]] <- cancerListPatients
pbmc[["Adjacent_Tumor"]] <- cancerListPatients2
pbmc[["residual"]]=t(abs(df.residual))
names(cancerListPatients) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, cancerListPatients)

a=df.residual*df.residual
pbmc[["residualsquared"]]=t(a)
pbmc[["residual"]]=abs(df.residual[,1])
pbmc[["pearson_residual"]]=abs(df.residual[,2])
pbmc[["studentized_residual"]]=abs(df.residual[,3])
pbmc[["standardized_residual"]]=abs(df.residual[,4])
pbmc[["exclusion"]]=exclusionlabel
pbmc[["exclusion_rsq"]]=exclusionlabelrsq1
pbmc[["FOXM1fit_forcedfit"]]=FOXM1fitrep
pbmc[["FOXM1fit"]]=FOXM1fitrep_exp
pbmc[["overall_survival"]]=surv_all[,1]
pbmc[["normalized_survival"]]=surv_all[,2]
pbmc[["median_survival"]]=surv_median_label
pbmc[["median_survival_zeros"]]=surv_median_zeros_label

pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
ElbowPlot(pbmc)
DimHeatmap(pbmc, dims = 13:18, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, reduction="umap")

VizDimLoadings(pbmc,reduction = "pca",dims=1:4)

DimPlot(pbmc, reduction = "umap")
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:20)

DimPlot(pbmc, reduction = "umap",group.by = "cancer",label = TRUE,repel = TRUE)
DimPlot(pbmc, reduction = "umap",group.by = "cancer",label=TRUE,repel = TRUE)      

DimPlot(pbmc, reduction = "pca",dims=c(3,4),group.by = "cancer",label=TRUE,repel = TRUE)

DimPlot(pbmc, reduction = "umap",group.by = "cancer",label=TRUE,repel = TRUE)

DimPlot(pbmc, reduction = "umap",group.by = "cancer",label=TRUE,repel = TRUE,label.size = 6)
DimPlot(pbmc, reduction = "umap",group.by = "cancer",label=TRUE,repel = TRUE,label.size = 6)+NoLegend()
DimPlot(pbmc, reduction = "umap",group.by = "Adjacent_Tumor")


DimPlot(pbmc, reduction = "umap",group.by = "exclusion",label = TRUE)
FeaturePlot(pbmc, features = "pearson_residual",order = TRUE)
FeaturePlot(pbmc, features = c("FOXM1","TOP2A","LMNB1","LMNA","COL4A1","COL4A2","ACTA2","COL1A1","COL1A2"))
plot1=DimPlot(pbmc, reduction = "umap",group.by = "cancer")+NoLegend()
plot2=DimPlot(pbmc, reduction = "umap",group.by = "exclusion")
plot3=FeaturePlot(pbmc, features = "residual",cols=c("red","green"),order = TRUE)
plot4=FeaturePlot(pbmc, features = "pearson_residual",order = TRUE)
plot5=FeaturePlot(pbmc, features = "studentized_residual",order = TRUE)
plot6=FeaturePlot(pbmc, features = "standardized_residual",order = TRUE)
plot1+plot2+plot3+plot4+plot5+plot6
plot1+plot2+plot3+plot7
plot1+plot2+plot7+plot14
plot1+plot2
plot1+plot8+plot3+plot7
plot1+plot2+plot7+plot11
plot7=FeaturePlot(pbmc, reduction="umap",features = "FOXM1fit",cols=c("green","red"))
plot11=FeaturePlot(pbmc, reduction="umap",features = "FOXM1fit_exp",order=TRUE)
plot8=DimPlot(pbmc, reduction = "umap",group.by = "exclusion_rsq",cols=c("lightseagreen","indianred1"),label = TRUE,repel = TRUE)
plot9=FeaturePlot(pbmc, reduction="umap",features = "overall_survival",order=TRUE)
plot10=FeaturePlot(pbmc, reduction="umap",features = "normalized_survival",order = TRUE)
plot11=FeaturePlot(pbmc, reduction="umap",features = "median_survival",order=TRUE,cols=c("blue","red","red4"),label = TRUE,repel = TRUE)
plot16=FeaturePlot(pbmc, reduction="umap",features = "median_survival",order=TRUE,cols=c("blue","red"),label = TRUE,repel = TRUE)

plot15=FeaturePlot(pbmc, reduction="umap",features = "median_survival",order=TRUE,cols=wes_palette("GrandBudapest1",n=3),label = TRUE,repel = TRUE)
plot15=FeaturePlot(pbmc, reduction="umap",features = "median_survival",order=TRUE,cols=c("rosybrown1","indianred1","red1"),label = TRUE,repel = TRUE)

plot14=FeaturePlot(pbmc, reduction="umap",features = "median_survival_zeros",order = TRUE,cols=c("green","red"))
plot12=FeaturePlot(pbmc, features = "residualsquared",cols=c("red","green"),order=TRUE)
DimPlot(pbmc, reduction = "umap",group.by = "exclusion_rsq")
plot1+plot8+plot14+plot11
plot1+plot11+plot15+plot16

plot8+plot11
plot1+plot7+plot9+plot10
pbmc$CellType <- Idents(pbmc)
Idents(pbmc) <- "cancer"
c0 <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25,logfc.threshold = 0.25)
c1 <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25,logfc.threshold = 0.25)
c2 <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25,logfc.threshold = 0.25)
c3 <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25,logfc.threshold = 0.25)
c4 <- FindMarkers(pbmc, ident.1 = 4, min.pct = 0.25,logfc.threshold = 0.25)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
top10=c("COL1A1","COL1A2","ACTA2","COL4A1","MKI67","COL4A2","LMNB1","FOXM1")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


write.table(pbmc[["pca"]], file = "PCA_BULK_ALL", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
