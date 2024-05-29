#
## Loading packages
#

library(tidyverse)
library(DESeq2)
library(survminer)
library(survival)
library(nlme)
library(ggalluvial)
library(ggrepel)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(forestplot)
library(biomaRt)
library(fst)
library(stringr)
library(org.Hs.eg.db)
library(DescTools)
library(scales)



#
##  Calculation of ORACLE risk-score
#

#Filter and normalised count data
#gene expression at least 1 TPM in at least 20% of samples, variance stabilising tranformation

TRACERx_count <- read.csv("2022-11-06_tx421_count.csv",row.names = 1)
TRACERx_TPM <- read.csv("2022-11-06_tx421_TPM.csv", row.names = 1)

#
n_TPM <- 1
tmp <- TRACERx_TPM > n_TPM   # at least n TPM
tumour_sample_n <- ncol(TRACERx_TPM)
tpm_threshold <- 0.2*tumour_sample_n   #in at least n% of samples

#keep genes
keep <- rowSums(tmp) >= tpm_threshold
TRACERx_filt_count <- TRACERx_count[keep,]

#prepare DESeq Dataset
coldata<-data.frame(sample = colnames(TRACERx_filt_count),row.names = colnames(TRACERx_filt_count))
TRACERx_filt_count <- round(TRACERx_filt_count)
dds <- DESeqDataSetFromMatrix(countData = TRACERx_filt_count,
                              colData = coldata,
                              design = ~1)

#vst normalisation
vsd <- vst(dds, blind=TRUE, fitType="mean")
vsd <- assay(vsd)

#calculate ORACLE risk score
Supp_Table_5_ORACLE <- read.csv("Supp_Table_5_ORACLE.csv")
Tx421_full_riskscore <- data.frame(sample_name = colnames(vsd) , ORACLE_riskscore = colSums(vsd[Supp_Table_5_ORACLE$Gene.Symbol,]*Supp_Table_5_ORACLE$Model.Coefficient) )


#
##  Batch correction
#

#Perform batch correction for risk scores
cutoff <- read.csv("2018-11-07_oracle_de.novo_cut.off.csv")
cutoff <- cutoff$RiskScore[1]

#Tx100 risk scores
tracerx_oracle_rs <- read.csv("data_20230810\\tracerx_oracle_rs_md.csv")
tracerx_oracle_rs$sample_name <- gsub("-",tracerx_oracle_rs$sample_name, replacement = ".")

#compare between Tx100 and Tx421
compare.df <- left_join(Tx421_full_riskscore,tracerx_oracle_rs[,c("sample_name","RiskScore")])
compare.df <- compare.df[which(!is.na(compare.df$RiskScore)),]
colnames(compare.df)[2:3] <- c("Tx421","Tx100")


#plot non-corrected riskscore with Tx100
md <- lm(data = compare.df, Tx100~Tx421)
r <- summary(md)
ggplot(compare.df,aes(Tx421,Tx100)) + geom_point() + geom_smooth(method = "lm",se=F) + 
  ylab("2019 Risk Score") + xlab("2023 Risk Score") + ggtitle(paste0("n_TPM=",n_TPM, ", n_sample=",85)) +
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  scale_x_continuous(expand = c(0,0),limits = c(8,12)) +
  geom_hline(yintercept = cutoff,lty="dashed") + geom_vline(xintercept = cutoff,lty="dashed") +
  geom_text(aes(x=8.8, y=11.5, label = paste0("R-sq=",signif(r$r.squared,digits = 3), "\np=",signif(r$coefficients[2,4],digits = 2)))) +
  geom_text(aes(x=8.8, y=11.1, label = paste0("Y=",signif(r$coefficient[2,1],digits = 3), "X",signif(r$coefficients[1,1],digits = 3)))) +
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#corrected riskscore
Tx421_full_riskscore$ORACLE_riskscore <- Tx421_full_riskscore$ORACLE_riskscore * r$coefficients[2,1] + r$coefficients[1,1]


#
##  RNA signature risk-scores
#

#Obtain count data on full list of genes
TRACERx_count <- read.csv("2022-11-06_tx421_count.csv",row.names = 1)

#prepare DESeq Dataset
coldata<-data.frame(sample = colnames(TRACERx_count),row.names = colnames(TRACERx_count))
TRACERx_count <- round(TRACERx_count)
dds <- DESeqDataSetFromMatrix(countData = TRACERx_count,
                              colData = coldata,
                              design = ~1)

#vst normalisation
vsd <- vst(dds, blind=TRUE, fitType="mean")
vsd <- assay(vsd)


#import signatures
Li_JAMA <- read.csv("Li_JAMA_2017.csv")
Wang <- read.csv("Wang_Front_Immuno_2022.csv")
Zhao<-read.csv("Zhao_LC_2020.csv")
Song <- read.csv("Song_Sci Rep_2022.csv")
Jin <- read.csv("Jin_J Immunol Res_2022.csv")
LiFeng <- read.csv("FengLi_Sci Rep_2022.csv")

#Li et al. JAMA Oncol 2017
mat <- vsd[which( rownames(vsd) %in% c(Li_JAMA$IRG1,Li_JAMA$IRG2) ), ] %>% as.data.frame()
Li_JAMA$gene_pair <- paste(Li_JAMA$IRG1,Li_JAMA$IRG2,sep="_")

IRGP_score <- matrix(nrow = length(Li_JAMA$gene_pair), ncol = ncol(mat),
                     dimnames = list(Li_JAMA$gene_pair, colnames(mat))) %>% as.data.frame
IRG1 <- NA
IRG2 <- NA
for(i in 1:ncol(IRGP_score)){
  for(k in 1:nrow(Li_JAMA)){
    IRG1 <- mat[which(rownames(mat) == Li_JAMA$IRG1[k]),i];
    IRG2 <- mat[which(rownames(mat) == Li_JAMA$IRG2[k]),i];
    IRGP_score[k,i] <- ifelse(IRG1 < IRG2,1,0)
  }
}

Li_riskscore <- data.frame(Li_RS = colSums(IRGP_score[Li_JAMA$gene_pair,] * Li_JAMA$Coef), sample_name = colnames(IRGP_score))


#Wang et al. Front Immunol 2022
mat <- vsd[which( rownames(vsd) %in% Wang$Gene ), ] %>% as.data.frame()
Wang_riskscore <- data.frame(Wang_RS = colSums(mat[Wang$Gene,] * Wang$Coefficient), sample_name = colnames(mat))


#Zhao et al. Lung cancer 2020
mat <- vsd[which( rownames(vsd) %in% Zhao$Gene ), ] %>% as.data.frame()
Zhao_riskscore <- data.frame(Zhao_RS = colSums(mat[Zhao$Gene,] * Zhao$Coefficient), sample_name = colnames(mat))


#Song et al. Sci Rep 2022
mat <- vsd[which( rownames(vsd) %in% Song$Genes ), ] %>% as.data.frame()
Song_riskscore <- data.frame(Song_RS = colSums(mat[Song$Genes,] * Song$coef), sample_name = colnames(mat))


#Jin et al. J Immunol Res 2022
mat <- vsd[which( rownames(vsd) %in% Jin$Gene ), ] %>% as.data.frame()
Jin_riskscore <- data.frame(Jin_RS = colSums(mat[Jin$Gene,] * Jin$beta), sample_name = colnames(mat))


#Feng Li et al. Sci Rep 2022
mat <- vsd[which( rownames(vsd) %in% LiFeng$Gene ), ] %>% as.data.frame()
LiFeng_riskscore <- data.frame(LiFeng_RS = colSums(mat[LiFeng$Gene,] * LiFeng$beta), sample_name = colnames(mat))


#
Tx421_full_riskscore <- left_join(Tx421_full_riskscore, Li_riskscore, by = "sample_name")
Tx421_full_riskscore <- left_join(Tx421_full_riskscore, Wang_riskscore, by = "sample_name")
Tx421_full_riskscore <- left_join(Tx421_full_riskscore, Zhao_riskscore, by = "sample_name")
Tx421_full_riskscore <- left_join(Tx421_full_riskscore, Song_riskscore, by = "sample_name")
Tx421_full_riskscore <- left_join(Tx421_full_riskscore, Jin_riskscore, by = "sample_name")
Tx421_full_riskscore <- left_join(Tx421_full_riskscore, LiFeng_riskscore, by = "sample_name")
Tx421_full_riskscore$sample_name <- gsub("[.]",Tx421_full_riskscore$sample_name, replacement = "-")

#filter for LUAD cohort and join with anonymous id - CRUK ID
Tx421_clinicopatho <- read.csv("Tx421_merged_clinicopathological_data_20220726.nooutcome.csv")
Tx421_LUAD <- Tx421_clinicopatho[which(Tx421_clinicopatho$Histology_per_region == "Invasive adenocarcinoma"),]
Tx421_LUAD$sample_name[which(Tx421_LUAD$sample_name == "LTX0615_SU_T1-R5")] <- "LTX0615_SU_T2-R5" #This is a mislabelled tumour which belong to T2 cluster
Tx421_LUAD <- Tx421_LUAD[which(grepl("SU_T1", Tx421_LUAD$sample_name)),] #select primary tumour

Tx421_riskscore <- left_join(Tx421_LUAD[,c("sample_name","cruk_id","sample_name_cruk","tumour_id_mphase_cruk")], Tx421_full_riskscore, by = "sample_name")
Tx421_riskscore <- Tx421_riskscore[which(!is.na(Tx421_riskscore$ORACLE_riskscore)),]


#
##  isolate TRACERx validation cohort
#

#separate region info
Tx421_riskscore$RegionalID <- apply(Tx421_riskscore, 1, FUN = function(x){ strsplit(x[1], split = "-")[[1]][2] })

#assign ORACLE risk-class
Tx421_riskscore$bin <- ifelse(Tx421_riskscore$ORACLE_riskscore > cutoff,"High","Low")
risk_class <- table(Tx421_riskscore$cruk_id,Tx421_riskscore$bin)
risk_class <- data.frame(High=as.matrix(risk_class)[,"High"], Low=as.matrix(risk_class)[,"Low"])
risk_class$ORACLE_class <- NA
risk_class$ORACLE_class <- ifelse(risk_class$High > 0, paste(risk_class$ORACLE_class, "High", sep=""), risk_class$ORACLE_class)
risk_class$ORACLE_class <- ifelse(risk_class$Low > 0, paste(risk_class$ORACLE_class, "Low", sep=""), risk_class$ORACLE_class)
risk_class$ORACLE_class <- gsub(x=risk_class$ORACLE_class, pattern="NA", replacement="")
risk_class$ORACLE_class <- gsub(x=risk_class$ORACLE_class, pattern="HighLow", replacement="Discordant")
risk_class$cruk_id <- rownames(risk_class)
Tx421_riskscore <- left_join(Tx421_riskscore,risk_class,by = "cruk_id")

#prepare TRACERx321 validation cohort - excluding patients used in Biswas et al. 2019
tracerx_oracle_rs <- read.csv("C:\\Users\\cola2318\\Dropbox\\UCL RA data\\Tx421\\github\\tracerx_oracle_rs.csv")
Tx100 <- unique(tracerx_oracle_rs$PublicationID)
Tx321_riskscore <- Tx421_riskscore[which(!Tx421_riskscore$cruk_id %in% Tx100),]



#
##  load data
#

### Figure 1

#Tx validation cohort risk-scores
Tx321_riskscore <- read.csv("20230727_TRACERx321_riskscores.csv")

#2019 cutoff
cutoff <- read.csv("2018-11-07_oracle_de.novo_cut.off.csv")
cutoff <- cutoff$RiskScore[1]

#load four metrics result
TSB_metric1 <- read.csv("Discordant_perc.csv")
TSB_metric2 <- read.csv("Cluster_concordance_AUC.csv")
TSB_metric3 <- read.csv("ExpVar_byHouseham.csv")
TSB_metric4 <- read.csv("LeastBiopsy_TSB.csv")




### Figure 2

#Tx risk-score
Tx321_riskscore <- read.csv("20230727_TRACERx321_riskscores.csv")

#Tx clinical data
TRACERx_clinical <- readRDS("all_patient_df_20221109.RDS")
all_tumour_df_20220209 <- readRDS("all_tumour_df_20220209.RDS")
all_tumour_df_20220209$tumour_id_mphase_cruk <- apply(all_tumour_df_20220209, 1 , FUN = function(x){ if(grepl("Cluster", x[12])) {x[12] <- gsub(".*-", x[12], replacement = paste0(x[4] , "_"))} else {x[12] <- x[4]} } )

#overall survival data
Tx321_Survdata <- left_join(Tx321_riskscore, TRACERx_clinical[,c("cruk_id","cens_os","os_time")], by = "cruk_id")

#join clinicopathological
Tx321_Survdata <- left_join(Tx321_Survdata, all_tumour_df_20220209[,c("tumour_id_mphase_cruk","age","sex","pack_years_calculated","pTNMStage_v8","adjuvant_treatment_YN","IASLC_LUAD_grade")], by = "tumour_id_mphase_cruk")
Tx321_Survdata$pTNMStage_v8[which(grepl("CRUK0704",Tx321_Survdata$tumour_id_mphase_cruk))] <- "2b" #collision tumour

#censor survival at 5 years
Tx321_Survdata$os_time <- Tx321_Survdata$os_time/365
Tx321_Survdata$cens_os[which(Tx321_Survdata$os_time >5)]<-0
Tx321_Survdata$os_time[which(Tx321_Survdata$os_time>5)] <- 5

#survival data for each patient
Tx321_tumour_Survdata <- Tx321_Survdata[which(!duplicated(Tx321_Survdata$cruk_id)),]
Tx321_tumour_Survdata$ORACLE_class <- factor(Tx321_tumour_Survdata$ORACLE_class, levels = c("Low","Discordant","High"))

#calculate ORACLE mean score
tmp <- aggregate(Tx321_Survdata$ORACLE_riskscore, by = list(cruk_id = Tx321_Survdata$cruk_id), mean)
colnames(tmp)[2] <- "ORACLE_mean"
Tx321_tumour_Survdata <- left_join(Tx321_tumour_Survdata, tmp, by = "cruk_id")

#factor levels
Tx321_tumour_Survdata$sex <- factor(Tx321_tumour_Survdata$sex , levels = c("Male","Female"))
Tx321_tumour_Survdata$adjuvant_treatment_YN <- factor(Tx321_tumour_Survdata$adjuvant_treatment_YN , levels = c("No adjuvant","Adjuvant"))
Tx321_tumour_Survdata$TNM_combine <- NA
Tx321_tumour_Survdata$TNM_combine[which(grepl("1",Tx321_tumour_Survdata$pTNMStage_v8 ))] <- "I"
Tx321_tumour_Survdata$TNM_combine[which(grepl("2",Tx321_tumour_Survdata$pTNMStage_v8 ))] <- "II"
Tx321_tumour_Survdata$TNM_combine[which(grepl("3",Tx321_tumour_Survdata$pTNMStage_v8 ))] <- "III"
Tx321_tumour_Survdata$TNM_combine <- factor(Tx321_tumour_Survdata$TNM_combine, levels = c("I","II","III"))
Tx321_tumour_Survdata$IASLC_LUAD_grade[which(Tx321_tumour_Survdata$IASLC_LUAD_grade == "Grade 1")] <- NA
Tx321_tumour_Survdata$IASLC_LUAD_grade <- factor(Tx321_tumour_Survdata$IASLC_LUAD_grade, levels = c("IMA","Grade 2","Grade 3"))


## 2C

#function for simulation - psuedo single biopsy cohort
create_surv <- function(signature){
  
  #
  RiskScore.df <- dplyr::select(sample.df, cruk_id, paste(signature,"RS",sep = "_"))
  
  #sample 1 region per patient tumour
  bootstrap.df <- aggregate(RiskScore.df[,paste(signature,"RS",sep = "_")] , by=list(cruk_id = RiskScore.df$cruk_id),FUN=function(x){sample(x,size=1,replace = T)})
  colnames(bootstrap.df)[2] <- "Boot_RS"
  
  #join os_time, cens_os
  surv.df <- left_join(sample.df,bootstrap.df,by="cruk_id")
  
  #use original riskscore for patient tumour with only 1 region
  surv.df[which(surv.df$cruk_id %in% names(which(table(surv.df$cruk_id)==1))),"Boot_RS"] <- surv.df[which(surv.df$cruk_id %in% names(which(table(surv.df$cruk_id)==1))),paste(signature,"RS",sep = "_")]
  
  #
  surv.df <- surv.df[which(!duplicated(surv.df$cruk_id)),c("cruk_id","Boot_RS")]
  
  #TRACERx risk ROC cut-off
  if(signature == "ORACLE"){
    cutoff <- read.csv("2018-11-07_oracle_de.novo_cut.off.csv")
    cutoff <- cutoff$RiskScore[1]
  }
  else{cutoff <- median(sample.df[,paste(signature,"RS",sep = "_")])}
  
  #bootstrap data frame
  surv.df$bin <- ifelse(surv.df[,"Boot_RS"] > cutoff, "High","Low")
  surv.df <- left_join(surv.df,Tx321_tumour_Survdata[,c("cruk_id","os_time","cens_os")],by="cruk_id")
  
  return(surv.df)
  
}

#function for survival association analysis - log-rank p, cox regression
bootsrap_surv <- function(survdf){
  
  survdiff <- survdiff(Surv(survdf$os_time,survdf$cens_os)~survdf$bin)
  logrank_p <- pchisq(survdiff$chisq,df=1,lower.tail = FALSE)
  logrank_p <- signif(logrank_p,digits = 2)
  
  #cox regression on bootstrap riskscore
  cox <- coxph(Surv(os_time, cens_os)~Boot_RS, data = survdf)
  cox_result <- summary(cox)
  cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])
  
  result.df <- data.frame(log_rank_p = logrank_p,HR = cox_result$HR, cox_pval = cox_result$P, cox_lci = cox_result$lower_ci, cox_uci = cox_result$upper_ci)
  
  return(result.df)
  
}




### Figure 3

## 3A-B

#lung-cancer-specific survival
Tx321_tumour_Survdata <- left_join(Tx321_tumour_Survdata, TRACERx_clinical[,c("cruk_id","cens_lung_specific","lung_specific_time")], by = "cruk_id")

#censor LCSS at 5 years
Tx321_tumour_Survdata$lung_specific_time <- Tx321_tumour_Survdata$lung_specific_time/365
Tx321_tumour_Survdata$cens_lung_specific[which(Tx321_tumour_Survdata$lung_specific_time >5)]<-0
Tx321_tumour_Survdata$lung_specific_time[which(Tx321_tumour_Survdata$lung_specific_time>5)] <- 5


## 3C

#Mascaux dataset
Mascaux <- read.delim("GSE33479.txt.gz", as.is=T, check.names=FALSE)
gene.annot <- read.delim("GPL6480-9577.txt", as.is=T, check.names=FALSE,skip=17)
patient.annot <- read.table("GSE33479_series_matrix.txt.gz", header = TRUE,fill = TRUE,skip=25)

#ORACLE genes
Supp_Table_5_ORACLE <- read.csv("Supp_Table_5_ORACLE.csv")

#remove non-matched probes
Mascaux <- Mascaux[,-c(124:127)]
Mascaux <- Mascaux[-which(is.na(Mascaux$ID_REF)),]
colnames(Mascaux)[1] <- "ID"

#match gene names to probes
Mascaux <- left_join(Mascaux,gene.annot[,c("ID","GENE_SYMBOL")], by="ID")
Mascaux <- Mascaux[-which(Mascaux$GENE_SYMBOL == ""),]

#paste patient ID
colnames(patient.annot) <- paste(patient.annot[15,], colnames(patient.annot),sep = "-")

#duplicated genes - max value
Mascaux_ORACLE <- Mascaux[which(Mascaux$GENE_SYMBOL %in% Supp_Table_5_ORACLE$Gene.Symbol),]
Mascaux_ORACLE <- aggregate(Mascaux_ORACLE[,-which(colnames(Mascaux_ORACLE) == "GENE_SYMBOL")],by=list(Gene = Mascaux_ORACLE$GENE_SYMBOL),max)
Mascaux_ORACLE <- column_to_rownames(Mascaux_ORACLE, var = "Gene")
Mascaux_ORACLE <- Mascaux_ORACLE[,-1]


## 3D
#load ORACLE clinical data - retreive recurrence samples
Tx421_riskscore <- read.csv("20230727_TRACERx421_riskscores.csv")
Tx421_clinicopatho <- read.csv("Tx421_merged_clinicopathological_data_20220726.nooutcome.csv")
TRACERx_clinical <- readRDS("all_patient_df_20221109.RDS")

#primary-metastasis phylogenies
seed.region <- read.table("seedingRegionInfo.txt",header = T,sep = "\t")


## 3E-F

#disease-free survival
Relapse_survdata <- dplyr::select(Tx321_riskscore, ORACLE_class, ORACLE_riskscore,cruk_id,tumour_id_mphase_cruk)
Relapse_survdata <- Relapse_survdata[which(!duplicated(Relapse_survdata$cruk_id)),]
Relapse_survdata <- left_join(Relapse_survdata, TRACERx_clinical[,c("cruk_id","dfs_time","cens_dfs","dfs_time_any_event","Relapse_cat_new","first_dfs_any_event")], by = "cruk_id")
Relapse_survdata <- left_join(Relapse_survdata, all_tumour_df_20220209[,c("tumour_id_mphase_cruk","age","sex","pack_years_calculated","pTNMStage_v8","adjuvant_treatment_YN")], by = "tumour_id_mphase_cruk")

# For DFS, censor 4 patients who are currently marked as "recurrence" but uncertain whether the recurrence is from 1st primary or 2nd primary (if from 2nd primary, these cases should NOT be marked as recurrence in TRACERx protocol)
#  "CRUK0512", "CRUK0373","CRUK0428","CRUK0511" : currently marked as recurrence, after 2nd primary cancer was confirmed
Relapse_survdata$cens_dfs[which(Relapse_survdata$cruk_id %in% c("CRUK0512","CRUK0373","CRUK0428","CRUK0511"))] <- 0
Relapse_survdata$dfs_time[which(Relapse_survdata$cruk_id %in% c("CRUK0512","CRUK0373","CRUK0428","CRUK0511"))] <- Relapse_survdata$dfs_time_any_event[which(Relapse_survdata$cruk_id %in% c("CRUK0512","CRUK0373","CRUK0428","CRUK0511"))]

#censor DFS at 5 years
Relapse_survdata$dfs_time_any_event <- Relapse_survdata$dfs_time_any_event/365
Relapse_survdata$dfs_time <- Relapse_survdata$dfs_time/365
Relapse_survdata$cens_dfs[which(Relapse_survdata$dfs_time >5)]<-0
Relapse_survdata$dfs_time[which(Relapse_survdata$dfs_time>5)] <- 5

#calculate ORACLE mean risk-score
tmp <- aggregate(Relapse_survdata$ORACLE_riskscore, by = list(cruk_id = Relapse_survdata$cruk_id), mean)
colnames(tmp)[2] <- "ORACLE_mean"
Relapse_survdata <- left_join(Relapse_survdata, tmp, by = "cruk_id")



### Figure 4

#CCLE test result
IC50_results <- read.csv("CCLE_result_table_20230308.csv",row.names = 1)



### Figure 5

#TRACERx clinicopathological risk factors + genetic evolutionary metrics
Features.df <- read.csv("Source Data ED8A.csv",check.names = F)

#TRACERx evolutionary biomarkers
Tx_biomarkers <- read.csv("Source Data Fig5B.csv")

#function for glm test
biomarker.glm <- function(FU_time){
  
  tmp <- dplyr::select(Tx_biomarker_avail, cruk_id,SCNA_ITH,STAS,ctDNA_status,subclonal_wgd, RecentSubclExpansionScore ,ORACLE_scaled,cens_os,os_time)
  
  #censor any patient exceed FU time
  tmp$cens_os[which(tmp$os_time >FU_time)]<-0
  tmp$os_time[which(tmp$os_time>FU_time)] <- FU_time
  
  #removed censored patients within the FU time
  if(FU_time > 1){
    tmp <- tmp[-which(tmp$os_time < FU_time & tmp$cens_os==0),]
  }
  
  #Create glm result data frame for all biomarker combinations
  colnames(tmp)[2:7] <- c("SCNA_ITH","STAS","ctDNA","S_WGD","RSE","ORACLE")
  
  glm.result.df <- data.frame(variable = c("SCNA_ITH","STAS","ctDNA","S_WGD","RSE","ORACLE"),Explained_var = NA)
  
  
  #
  for(i in 1:nrow(glm.result.df)){
    mod <- glm(data=tmp,as.formula(paste0("cens_os ~",glm.result.df$variable[i])) )
    glm.result.df$Explained_var[i] <- PseudoR2(mod)*100
  }
  
  glm_time <- glm.result.df
  glm_time$OS_time <- FU_time
  glm_time$event <- paste0(glm_time$OS_time,"(" ,table(tmp$cens_os)["1"],")")
  
  return(glm_time)
  
} 



### Extended Data Fig 2
all_tumour_df_20220209 <- readRDS("all_tumour_df_20220209.RDS")
all_tumour_df_20220209$tumour_id_mphase_cruk <- apply(all_tumour_df_20220209, 1 , FUN = function(x){ if(grepl("Cluster", x[12])) {x[12] <- gsub(".*-", x[12], replacement = paste0(x[4] , "_"))} else {x[12] <- x[4]} } )



### Extended Data Fig 3

#load signature genes
Li_JAMA <- read.csv("Li_JAMA_2017.csv")
Wang <- read.csv("Wang_Front_Immuno_2022.csv")
Zhao<-read.csv("Zhao_LC_2020.csv")
Song <- read.csv("Song_Sci Rep_2022.csv")
Jin <- read.csv("Jin_J Immunol Res_2022.csv")
LiFeng <- read.csv("FengLi_Sci Rep_2022.csv")

#load LUAD sample expression matrix
LUAD_vsd <- read.csv("tracerx_primaryLUAD_vsd.csv", row.names = 1)

#prepare expression matrix for each signature
mat1 <- LUAD_vsd[Supp_Table_5_ORACLE$Gene.Symbol,] #ORACLE
mat2 <- LUAD_vsd[which(rownames(LUAD_vsd) %in% c(Li_JAMA$IRG1,Li_JAMA$IRG2)), ] #Li
mat3 <- LUAD_vsd[which(rownames(LUAD_vsd) %in% Wang$Gene), ] #Wang
mat4 <- LUAD_vsd[which(rownames(LUAD_vsd) %in% Zhao$Gene), ] #Zhao
mat5 <- LUAD_vsd[which(rownames(LUAD_vsd) %in% Song$Genes), ] #Song
mat6 <- LUAD_vsd[which(rownames(LUAD_vsd) %in% Jin$Gene), ] #Jin
mat7 <- LUAD_vsd[which(rownames(LUAD_vsd) %in% LiFeng$Gene), ] #Feng Li

#function for calculating concordant rate in the same cluster
clusterSigs <- function(mat,k){
  
  #hclust
  tmp <- scale(t(mat),center = TRUE)
  clust <- hclust(dist(tmp,method = "manhattan"),method = "ward.D2")
  
  #cut into k clusters
  clust_data <- data.frame(cluster=cutree(clust,k), sample_name = names(cutree(clust,k)))
  
  #paste IDs to cluster
  clust_data$Shorter_ID <- gsub("_.*",rownames(clust_data),replacement = "")
  for(i in 1:nrow(clust_data)){
    clust_data$cluster[i] <- paste(clust_data$Shorter_ID[i],clust_data$cluster[i],sep = "-")
  }
  
  #assign concordant or discordant cluster
  clust_data$result <- NA
  for(i in 1:nrow(clust_data)){
    if(length(clust_data$cluster[which(clust_data$Shorter_ID==clust_data$Shorter_ID[i])]%>%unique()) == 1){
      clust_data$result[i] <- "Concordant"
    }
    else {clust_data$result[i] <- "Discordant"}
  }
  
  #calculate percentage
  clust_data_summary <- clust_data[-which(duplicated(clust_data$Shorter_ID)),]
  return(100*table(clust_data_summary$result)["Concordant"]/nrow(clust_data_summary))
}

#function to create df for samples fall in which cluster group when cut into k clusters
clusterdf <- function(mat, ht, k){
  
  tmp <- scale(t(mat),center = TRUE,scale = T)
  df <- data.frame(sample_name=NA,clust=NA)
  
  #from 1 to k cluster
  for(i in 1:k){
    a <- data.frame(sample_name = colnames(t(tmp))[column_order(draw(ht))[[i]]], clust=i)
    df <- rbind(df,a)
  }
  
  df <- df[-1,]
  df$cruk_id <- gsub("_.*",df$sample,replacement = "")
  rownames(df) <- df$sample
  df <- dplyr::select(df,clust)
  return(df)
}

#function to create heatmap
clust_heatmap <- function(mat){
  
  tmp <- scale(t(mat),center = TRUE,scale = T)
  
  ##top heatmap
  ht_gene <- Heatmap(t(tmp),
                     clustering_distance_columns =  dist(tmp,method = "manhattan"),
                     clustering_method_columns = "ward.D2",
                     cluster_rows = FALSE,
                     name = "ht1",
                     row_dend_width = unit(15, "mm"),
                     row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=2),
                     col = colorRamp2(seq(-3,3), viridis(7)),
                     show_row_names = TRUE,show_column_names = F,row_names_side = "left",
                     heatmap_legend_param = list(title = NULL, color_bar = "continuous")
  )
  
  ##cluster number annotation 
  ht2 <- Heatmap(t(tmp),
                 clustering_distance_columns =  dist(tmp,method = "manhattan"),
                 clustering_method_columns = "ward.D2",
                 cluster_rows = FALSE,
                 name = "ht2",
                 column_split = 2,
                 row_dend_width = unit(15, "mm"),
                 row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=2),
                 col = colorRamp2(seq(-3,3), viridis(7)),
                 show_row_names = TRUE,show_column_names = F,row_names_side = "left",
                 heatmap_legend_param = list(title = NULL, color_bar = "continuous")
  )
  
  df2 <- clusterdf(mat, ht2,2)
  df2 <- data.frame(cluster = df2[colnames(mat),])
  rownames(df2) <- colnames(mat)
  colors <- structure(colorRampPalette(brewer.pal(n=11, name="RdYlBu")[-c(1, 11)])(2), names = c(1:2))
  ha2 <-HeatmapAnnotation("2" = df2[colnames(mat), ],height = unit(3,"cm"),border = T,show_legend = F,
                          show_annotation_name = T,col = list("2"=colors),annotation_name_side="left")
  
  #
  ht10 <- Heatmap(t(tmp),
                  clustering_distance_columns =  dist(tmp,method = "manhattan"),
                  clustering_method_columns = "ward.D2",
                  cluster_rows = FALSE,
                  name = "ht10",
                  column_split = 10,
                  row_dend_width = unit(15, "mm"),
                  row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=2),
                  col = colorRamp2(seq(-3,3), viridis(7)),
                  show_row_names = TRUE,show_column_names = F,row_names_side = "left",
                  heatmap_legend_param = list(title = NULL, color_bar = "continuous")
  )
  
  df10 <- clusterdf(mat, ht10,10)
  df10 <- data.frame(cluster = df10[colnames(mat),])
  rownames(df10) <- colnames(mat)
  colors <- structure(colorRampPalette(brewer.pal(n=11, name="RdYlBu")[-c(1, 11)])(10), names = c(1:10))
  ha10 <-HeatmapAnnotation("10" = df10[colnames(mat), ],height = unit(3,"cm"),border = T,
                           show_legend = F,show_annotation_name = T,col = list("10"=colors),annotation_name_side="left")
  
  #
  ht60 <- Heatmap(t(tmp),
                  clustering_distance_columns =  dist(tmp,method = "manhattan"),
                  clustering_method_columns = "ward.D2",
                  cluster_rows = FALSE,
                  name = "ht60",
                  column_split = 60,
                  row_dend_width = unit(15, "mm"),
                  row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=2),
                  col = colorRamp2(seq(-3,3), viridis(7)),
                  show_row_names = TRUE,show_column_names = F,row_names_side = "left",
                  heatmap_legend_param = list(title = NULL, color_bar = "continuous")
  )
  
  df60 <- clusterdf(mat, ht60,60)
  df60 <- data.frame(cluster = df60[colnames(mat),])
  rownames(df60) <- colnames(mat)
  colors <- structure(colorRampPalette(brewer.pal(n=11, name="RdYlBu")[-c(1, 11)])(60), names = c(1:60))
  ha60 <-HeatmapAnnotation("60" = df60[colnames(mat), ],height = unit(3,"cm"),border = T,
                           show_legend = F,show_annotation_name = T,col = list("60"=colors),annotation_name_side="left")
  
  ##patient heatmap
  idx <- data.frame(sample_name = colnames(mat))
  idx <- left_join(idx, Tx321_riskscore[,c("sample_name","cruk_id")], by ="sample_name")
  idx$value <- "Y"
  
  patient_ht <- as.data.frame(spread(idx,sample_name,value = value))
  rownames(patient_ht) <- patient_ht$cruk_id
  order <- data.frame(sample_name = colnames(mat)[column_order(draw(ht1))])
  order <- left_join(order, Tx321_riskscore[,c("sample_name","cruk_id")], by ="sample_name")
  patient_ht <- patient_ht[unique(order$cruk_id),-1]
  
  ht_patient <- Heatmap(patient_ht,
                        cluster_columns = F,
                        cluster_rows = FALSE,
                        name = "ht3",
                        na_col = "white",
                        show_heatmap_legend = F,
                        row_names_gp = gpar(fontsize=5),
                        col = structure("black",names="Y"), row_order = unique(order$cruk_id),
                        show_row_names = TRUE,show_column_names = FALSE,row_names_side = "left",width = 1
  )
  
  #combine heatmap
  ht_list <- ht1%v%ha2%v%ha10%v%ha60%v%ht3
  ht_list
  print(decorate_heatmap_body("ht1", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))});
        decorate_heatmap_body("ht3", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))});
        decorate_heatmap_body("ht3", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))});
        decorate_heatmap_body("ht3", {for (i in 1:122){grid.lines(x=unit(c(1,0),"npc"),y=unit(c(i/122,i/122),"npc"), gp = gpar(lty="dotted", col="gray75", lwd=0.75))};grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})
  )
  
  
  ##pie chart for each cluster
  df2$sample_name <- rownames(df2)
  df2 <- left_join(df2, Tx321_riskscore[,c("sample_name","cruk_id")])
  df2$group <- apply(df2,1,FUN = function(x){paste(x[3],x[1],sep = "-")})
  
  #clust 2
  df2$result <- NA
  for(i in 1:nrow(df2)){
    if(length(unique(df2$group[which(df2$cruk_id==df2$cruk_id[i])])) == 1){
      df2$result[i] <- "Concordant"
    }
    else {df2$result[i] <- "Discordant"}
  }
  tmp <- df2[-which(duplicated(df2$cruk_id)),]
  pie.df <- data.frame(class = c("Concordant","Discordant"),perc=c(length(which(tmp$result=="Concordant"))/nrow(tmp),length(which(tmp$result=="Discordant"))/nrow(tmp)))
  ggplot(pie.df)+geom_col(aes(x=1,y=perc,fill=class))+coord_polar("y", start=0)+
    scale_fill_manual(values = c("black","azure4"))+
    theme_void() + theme(legend.position = "none")
  
  #clust 10
  df10$sample_name <- rownames(df10)
  df10 <- left_join(df10, Tx321_riskscore[,c("sample_name","cruk_id")])
  df10$group <- apply(df10,1,FUN = function(x){paste(x[3],x[1],sep = "-")})
  
  #
  df10$result <- NA
  for(i in 1:nrow(df10)){
    if(length(unique(df10$group[which(df10$cruk_id==df10$cruk_id[i])])) == 1){
      df10$result[i] <- "Concordant"
    }
    else {df10$result[i] <- "Discordant"}
  }
  tmp2 <- df10[-which(duplicated(df10$cruk_id)),]
  pie.df <- data.frame(class = c("Concordant","Discordant"),perc=c(length(which(tmp2$result=="Concordant"))/nrow(tmp2),length(which(tmp2$result=="Discordant"))/nrow(tmp2)))
  ggplot(pie.df)+geom_col(aes(x=1,y=perc,fill=class))+coord_polar("y", start=0)+
    scale_fill_manual(values = c("black","azure4"))+
    theme_void() + theme(legend.position = "none")
  
  #clust 60
  df60$sample_name <- rownames(df60)
  df60 <- left_join(df60, Tx321_riskscore[,c("sample_name","cruk_id")])
  df60$group <- apply(df60,1,FUN = function(x){paste(x[3],x[1],sep = "-")})
  
  #
  df60$result <- NA
  for(i in 1:nrow(df60)){
    if(length(unique(df60$group[which(df60$cruk_id==df60$cruk_id[i])])) == 1){
      df60$result[i] <- "Concordant"
    }
    else {df60$result[i] <- "Discordant"}
  }
  tmp3 <- df60[-which(duplicated(df60$cruk_id)),]
  pie.df <- data.frame(class = c("Concordant","Discordant"),perc=c(length(which(tmp3$result=="Concordant"))/nrow(tmp3),length(which(tmp3$result=="Discordant"))/nrow(tmp3)))
  ggplot(pie.df)+geom_col(aes(x=1,y=perc,fill=class))+coord_polar("y", start=0)+
    scale_fill_manual(values = c("black","azure4"))+
    theme_void() + theme(legend.position = "none")
  
  ##barplot
  rownames(tmp) <- tmp$cruk_id
  tmp <- tmp[rownames(patient_ht),]
  tmp$order <- c(nrow(tmp):1)
  ggplot(tmp)+geom_tile(aes(x=fct_reorder(cruk_id,order),y=1,fill=result))+coord_flip()+
    scale_fill_manual(values = c("black","azure4"))+
    theme_void() + theme(legend.position = "none")
  
  rownames(tmp2) <- tmp2$cruk_id
  tmp2 <- tmp2[rownames(patient_ht),]
  tmp2$order <- c(nrow(tmp2):1)
  ggplot(tmp2)+geom_tile(aes(x=fct_reorder(cruk_id,order),y=1,fill=result))+coord_flip()+
    scale_fill_manual(values = c("black","azure4"))+
    theme_void() + theme(legend.position = "none")
  
  rownames(tmp3) <- tmp3$cruk_id
  tmp3 <- tmp3[rownames(patient_ht),]
  tmp3$order <- c(nrow(tmp3):1)
  ggplot(tmp3)+geom_tile(aes(x=fct_reorder(cruk_id,order),y=1,fill=result))+coord_flip()+
    scale_fill_manual(values = c("black","azure4"))+
    theme_void() + theme(legend.position = "none")
  
}



### Extended Data Fig 4

## 4C-D

#function to calculate expression sd
ExpSD.df <- function(mat){
  
  df <- as.data.frame(t(mat))
  df$PatientID <- gsub("_.*",rownames(df),replacement = "")
  
  #calculate mean expression
  df_mean <- aggregate(df[,-ncol(df)],by=list(df$PatientID), FUN = function(x){ mean(x,na.rm=T)})
  
  #calculate sd of expression
  df_sd <- aggregate(df[,-ncol(df)],by=list(df$PatientID),FUN = function(x){ sd(x,na.rm=T)})
  
  #remove single-region samples
  df_sd <- df_sd[,which(!grepl("NA",colnames(df_sd)))]
  
  #calculate average sd
  rownames(df_sd) <- df_sd$Group.1
  df_sd <- df_sd[,-1]
  annot.df <- data.frame(mean_var = colSums(df_sd)/nrow(df_sd))
  
  #descending order
  annot.df <- arrange(annot.df,desc(mean_var))
  df_sd <- df_sd[,rownames(annot.df)]
  
  return(df_sd)
}

#function to calculate mean variation
MeanSD.df <- function(mat,sig){
  
  df <- as.data.frame(t(mat))
  df$PatientID <- gsub("_.*",rownames(df),replacement = "")
  
  #calculate mean expression
  df_mean <- aggregate(df[,-ncol(df)],by=list(df$PatientID), FUN = function(x){ mean(x,na.rm=T)})
  
  #calculate sd of expression
  df_sd <- aggregate(df[,-ncol(df)],by=list(df$PatientID),FUN = function(x){ sd(x,na.rm=T)})
  
  #calculate sd of expression
  df_sd <- df_sd[,which(!grepl("NA",colnames(df_sd)))]
  
  #remove single-region samples
  rownames(df_sd) <- df_sd$Group.1
  df_sd <- df_sd[,-1]
  
  #calculate average sd
  annot.df <- data.frame(mean_var = colSums(df_sd)/nrow(df_sd))
  annot.df <- arrange(annot.df,desc(mean_var))
  annot.df$signature <- sig
  
  return(annot.df)
}



### Extended Data Fig 5

#load microarray data
load("2019-04-17_microarray_expression.RData")
load("2019-04-17_microarray_survival.RData")



### Extended Data Fig 6

#GDSC IC50 data
LUAD_IC50 <- read.csv("C:\\Users\\cola2318\\Desktop\\UCL RA\\RA project\\CCLE\\LUAD_IC50_20220112.csv",row.names = 1)
df <- read_excel("C:\\Users\\cola2318\\Desktop\\UCL RA\\RA project\\CCLE\\LUAD_IC50_20220112.xlsx")
colnames(LUAD_IC50) <- colnames(df)[-1]



### Extended Data Fig 9

#load regional mutation data
Tx421_regionMut <- fst::read.fst("tx421_RegionMutTable.fst")

#load GISTIC results
HighLow_summary <- read.csv("20230620_concordant_oracle_highlow_cytoband_summary_wide.csv")

#filter for driver mutations
Tx421_regionMut <- Tx421_regionMut[which(Tx421_regionMut$PASS == "TRUE" & Tx421_regionMut$Is.present.region == "TRUE"),]
Tx421_regionMut <- Tx421_regionMut[which(Tx421_regionMut$DriverMut == "TRUE"),]

#fix names
Tx421_regionMut$RegionID <- gsub("_LN0",Tx421_regionMut$RegionID,replacement = "_LN")
Tx421_regionMut$RegionID <- gsub("_LN",Tx421_regionMut$RegionID,replacement = "_LN0")
Tx421_regionMut$RegionID <- gsub(":SU_T1.",Tx421_regionMut$RegionID,replacement = "_SU_T1-")
Tx421_regionMut$RegionID <- gsub(":SU_T2.",Tx421_regionMut$RegionID,replacement = "_SU_T2-")
Tx421_regionMut$RegionID <- gsub(":SU_LN",Tx421_regionMut$RegionID,replacement = "_SU_LN")



### Extended Data Fig 10

#load Uppsala LUSC
UPP_count_data <- read.delim("UPP LUSC\\GSE81089_readcounts_featurecounts.tsv.gz")
UPP_sample <- read.delim("UPP LUSC\\GSE81089_series_matrix.txt.gz",skip = 64)




#
## plot figures
#


### Figure 1

## 1B

#keep only multi-region samples
Tx321_samplingbias <- Tx321_riskscore[-which(Tx321_riskscore$High == 0 & Tx321_riskscore$Low == 1),]
Tx321_samplingbias <- Tx321_samplingbias[-which(Tx321_samplingbias$High == 1 & Tx321_samplingbias$Low == 0),]

#join stage info
Tx321_samplingbias$ORACLE_class <- factor(Tx321_samplingbias$ORACLE_class,levels = c("Low","Discordant","High"))

ggplot(Tx321_samplingbias,aes(x=fct_reorder(cruk_id, ORACLE_riskscore + as.numeric(ORACLE_class), .fun=mean), y = ORACLE_riskscore))+
  geom_point(aes(col = ORACLE_class),alpha = 0.5,pch=16,size = 5)+
  geom_hline(yintercept = cutoff,lty = "dotted")+
  geom_line(col = "black",size=0.5)+
  xlab("Patient ID")+
  ylab("ORACLE Risk score")+
  labs(color="Class")+
  scale_y_continuous(breaks = seq(8.5,11.5,1),expand = c(0,0),limits = c(8.5,11.7))+
  scale_x_discrete(expand = c(0.01,0.01)) +
  scale_color_manual(values = c("#3B4992FF","azure4","#EE0000FF"))+
  theme(legend.position = "none",panel.background = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))


## 1C

#dichotomise tumour region risk
Tx321_samplingbias$Li_bin <- ifelse(Tx321_samplingbias$Li_RS < median(Tx321_samplingbias$Li_RS), "Low","High")
Tx321_samplingbias$Wang_bin <- ifelse(Tx321_samplingbias$Wang_RS < median(Tx321_samplingbias$Wang_RS), "Low","High")
Tx321_samplingbias$Zhao_bin <- ifelse(Tx321_samplingbias$Zhao_RS < median(Tx321_samplingbias$Zhao_RS), "Low","High")
Tx321_samplingbias$Song_bin <- ifelse(Tx321_samplingbias$Song_RS < median(Tx321_samplingbias$Song_RS), "Low","High")
Tx321_samplingbias$Jin_bin <- ifelse(Tx321_samplingbias$Jin_RS < median(Tx321_samplingbias$Jin_RS), "Low","High")
Tx321_samplingbias$LiFeng_bin <- ifelse(Tx321_samplingbias$LiFeng_RS < median(Tx321_samplingbias$LiFeng_RS), "Low","High")

#assign patient risk-class for each signature
#Li JAMA Oncol 2017
tmp <- table(Tx321_samplingbias$cruk_id,Tx321_samplingbias$Li_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Li_class <- NA
tmp$Li_class <- ifelse(tmp$High > 0, paste(tmp$Li_class, "High", sep=""), tmp$Li_class)
tmp$Li_class <- ifelse(tmp$Low > 0, paste(tmp$Li_class, "Low", sep=""), tmp$Li_class)
tmp$Li_class <- gsub(x=tmp$Li_class, pattern="NA", replacement="")
tmp$Li_class <- gsub(x=tmp$Li_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Li_class <- factor(tmp$Li_class, levels = c("Low","Discordant","High"))
Tx321_samplingbias <- left_join(Tx321_samplingbias, tmp[,c("cruk_id","Li_class")], by = "cruk_id")

#Wang Front Immunol 2022
tmp <- table(Tx321_samplingbias$cruk_id,Tx321_samplingbias$Wang_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Wang_class <- NA
tmp$Wang_class <- ifelse(tmp$High > 0, paste(tmp$Wang_class, "High", sep=""), tmp$Wang_class)
tmp$Wang_class <- ifelse(tmp$Low > 0, paste(tmp$Wang_class, "Low", sep=""), tmp$Wang_class)
tmp$Wang_class <- gsub(x=tmp$Wang_class, pattern="NA", replacement="")
tmp$Wang_class <- gsub(x=tmp$Wang_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Wang_class <- factor(tmp$Wang_class, levels = c("Low","Discordant","High"))
Tx321_samplingbias <- left_join(Tx321_samplingbias, tmp[,c("cruk_id","Wang_class")], by = "cruk_id")

#Zhao Lung Cancer 2020
tmp <- table(Tx321_samplingbias$cruk_id,Tx321_samplingbias$Zhao_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Zhao_class <- NA
tmp$Zhao_class <- ifelse(tmp$High > 0, paste(tmp$Zhao_class, "High", sep=""), tmp$Zhao_class)
tmp$Zhao_class <- ifelse(tmp$Low > 0, paste(tmp$Zhao_class, "Low", sep=""), tmp$Zhao_class)
tmp$Zhao_class <- gsub(x=tmp$Zhao_class, pattern="NA", replacement="")
tmp$Zhao_class <- gsub(x=tmp$Zhao_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Zhao_class <- factor(tmp$Zhao_class, levels = c("Low","Discordant","High"))
Tx321_samplingbias <- left_join(Tx321_samplingbias, tmp[,c("cruk_id","Zhao_class")], by = "cruk_id")

#Song Sci Rep 2022
tmp <- table(Tx321_samplingbias$cruk_id,Tx321_samplingbias$Song_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Song_class <- NA
tmp$Song_class <- ifelse(tmp$High > 0, paste(tmp$Song_class, "High", sep=""), tmp$Song_class)
tmp$Song_class <- ifelse(tmp$Low > 0, paste(tmp$Song_class, "Low", sep=""), tmp$Song_class)
tmp$Song_class <- gsub(x=tmp$Song_class, pattern="NA", replacement="")
tmp$Song_class <- gsub(x=tmp$Song_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Song_class <- factor(tmp$Song_class, levels = c("Low","Discordant","High"))
Tx321_samplingbias <- left_join(Tx321_samplingbias, tmp[,c("cruk_id","Song_class")], by = "cruk_id")

#Jin J Immunol Res 2022
tmp <- table(Tx321_samplingbias$cruk_id,Tx321_samplingbias$Jin_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Jin_class <- NA
tmp$Jin_class <- ifelse(tmp$High > 0, paste(tmp$Jin_class, "High", sep=""), tmp$Jin_class)
tmp$Jin_class <- ifelse(tmp$Low > 0, paste(tmp$Jin_class, "Low", sep=""), tmp$Jin_class)
tmp$Jin_class <- gsub(x=tmp$Jin_class, pattern="NA", replacement="")
tmp$Jin_class <- gsub(x=tmp$Jin_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Jin_class <- factor(tmp$Jin_class, levels = c("Low","Discordant","High"))
Tx321_samplingbias <- left_join(Tx321_samplingbias, tmp[,c("cruk_id","Jin_class")], by = "cruk_id")

#Li Sci Rep 2022
tmp <- table(Tx321_samplingbias$cruk_id,Tx321_samplingbias$LiFeng_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$LiFeng_class <- NA
tmp$LiFeng_class <- ifelse(tmp$High > 0, paste(tmp$LiFeng_class, "High", sep=""), tmp$LiFeng_class)
tmp$LiFeng_class <- ifelse(tmp$Low > 0, paste(tmp$LiFeng_class, "Low", sep=""), tmp$LiFeng_class)
tmp$LiFeng_class <- gsub(x=tmp$LiFeng_class, pattern="NA", replacement="")
tmp$LiFeng_class <- gsub(x=tmp$LiFeng_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$LiFeng_class <- factor(tmp$LiFeng_class, levels = c("Low","Discordant","High"))
Tx321_samplingbias <- left_join(Tx321_samplingbias, tmp[,c("cruk_id","LiFeng_class")], by = "cruk_id")

#calculate risk-class frequency
TSB_summary <- Tx321_samplingbias[which(!duplicated(Tx321_samplingbias$cruk_id)),]
TSB <- rbind(table(TSB_summary$ORACLE_class)*100/nrow(TSB_summary),
             table(TSB_summary$Li_class)*100/nrow(TSB_summary),
             table(TSB_summary$Wang_class)*100/nrow(TSB_summary),
             table(TSB_summary$Zhao_class)*100/nrow(TSB_summary),
             table(TSB_summary$Song_class)*100/nrow(TSB_summary),
             table(TSB_summary$Jin_class)*100/nrow(TSB_summary),
             table(TSB_summary$LiFeng_class)*100/nrow(TSB_summary)) %>% as.data.frame()

#tidy up variable level
TSB$signature <- c("ORACLE","Li JAMA Oncol 2017","Wang Front Immunol 2022","Zhao Lung Cancer 2020","Song Sci Rep 2022","Jin J Immunol Res 2022","Li Sci Rep 2022")
TSB$signature <- factor(TSB$signature, levels = c("ORACLE" ,"Li JAMA Oncol 2017","Wang Front Immunol 2022","Zhao Lung Cancer 2020","Song Sci Rep 2022","Jin J Immunol Res 2022","Li Sci Rep 2022"))
TSB_melt <- melt(TSB)
TSB_melt$variable <- factor(TSB_melt$variable, levels = c("Low","High","Discordant"))

ggplot(TSB_melt, aes(x = "",y=value, fill=variable))+
  geom_col(width = 1, col = "black")+
  coord_polar("y", start=0)+
  ylab("")+ xlab("")+ 
  geom_text(aes(label = paste0(signif(value,digits = 2),"%")),position=position_stack(vjust=0.5),size=6)+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  facet_grid(~signature) +
  theme_void()+
  theme(legend.position = "none")


## 1D

#tidy up
TSB_metric1 <- dplyr::select(TSB_metric1, signature, Discordant)
TSB_metric2 <- dplyr::select(TSB_metric2, signature, AUC)
TSB_metric3 <- aggregate(TSB_metric3$mean_var, by = list(signature = TSB_metric3$signature), median) # median variation across signature genes

#ranking by each metric
TSB_metric1 <- arrange(TSB_metric1, Discordant)
TSB_metric1$rank <- 1:nrow(TSB_metric1)

TSB_metric2 <- arrange(TSB_metric2, desc(AUC))
TSB_metric2$rank <- 1:nrow(TSB_metric2)

TSB_metric3 <- arrange(TSB_metric3, x)
TSB_metric3$rank <- 1:nrow(TSB_metric3)

TSB_metric4 <- arrange(TSB_metric4, biop_need)
TSB_metric4$rank <- 1:nrow(TSB_metric4)

#unify column names
colnames(TSB_metric1)[2] <- "value"
colnames(TSB_metric2)[2] <- "value"
colnames(TSB_metric3)[2] <- "value"
colnames(TSB_metric4)[2] <- "value"

#round up
TSB_metric1$value <- paste0(signif(TSB_metric1$value,digits = 2), "%")
TSB_metric2$value <- signif(TSB_metric2$value,digits = 3)
TSB_metric3$value <- signif(TSB_metric3$value,digits = 2)
TSB_metric4$value <- signif(TSB_metric4$value,digits = 3)

#
TSB_metric1 <- dplyr::select(TSB_metric1,signature, rank,value)
TSB_metric2 <- dplyr::select(TSB_metric2,signature, rank,value)
TSB_metric3 <- dplyr::select(TSB_metric3,signature, rank,value)
TSB_metric4 <- dplyr::select(TSB_metric4,signature, rank,value)

#assign metric names
TSB_metric1$metric <- "Discordant rate"
TSB_metric2$metric <- "AUC"
TSB_metric3$metric <- "Expression SD"
TSB_metric4$metric <- "Least Biopsy"

All_TSB_metrics <- rbind(TSB_metric1, TSB_metric2, TSB_metric3, TSB_metric4)

#calculate mean rank
tmp <- aggregate(All_TSB_metrics$rank , by = list(signature = All_TSB_metrics$signature), mean)
tmp <- arrange(tmp, x)
tmp$rank <- 1:nrow(tmp)
colnames(tmp)[2] <- "value"
tmp$metric <- "mean"
tmp <- dplyr::select(tmp, signature, rank, value, metric)

#3 sigs with equal mean rank
tmp$rank[which(tmp$value == 4.5)] <- 3

#order by ranking
All_TSB_metrics <- rbind(All_TSB_metrics, tmp)
All_TSB_metrics$rank <- factor(All_TSB_metrics$rank,levels = c(7,6,5,4,3,2,1))
All_TSB_metrics$metric <- factor(All_TSB_metrics$metric,levels = c("Discordant rate","AUC","Expression SD","Least Biopsy","mean"))
All_TSB_metrics$signature <- factor(All_TSB_metrics$signature, levels = c("ORACLE","Li JAMA Oncol 2017","Wang Front Immunol 2022","Zhao Lung Cancer 2020","Song Sci Rep 2022","Jin J Immunol Res 2022","Li Sci Rep 2022"))
All_TSB_metrics$label <- ifelse(All_TSB_metrics$metric == "mean", "Yes","No")

ggplot(All_TSB_metrics, aes(metric, rank)) + 
  geom_point(aes(col=signature,size = label))+
  geom_line(aes(group = signature,col=signature),size=2) +
  xlab("") + ylab("Ranking") +
  geom_text(aes(x=metric, y=rank, label = value)) + 
  scale_size_manual(values = c(5,10)) +
  scale_color_manual(values = c("#A6CEE3","#1F78B4","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")) +
  theme(panel.background = element_blank(),panel.border = element_blank(),axis.ticks = element_blank())




### Figure 2

## 2A

j <- ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Tx321_tumour_Survdata),
                legend.labs = c( "Low","Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                legend="none", xlab="Time (years)", ylab="Overall Survival",
                pval = T,pval.size = 6,pval.coord = c(0, 0.05),
                ggtheme = theme(title = element_text(size = 18),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j


## 2B

MVA.cox.res <- coxph(formula = Surv(Tx321_tumour_Survdata$os_time, Tx321_tumour_Survdata$cens_os) ~ 
                       Tx321_tumour_Survdata$ORACLE_mean + 
                       Tx321_tumour_Survdata$sex +
                       Tx321_tumour_Survdata$age + 
                       Tx321_tumour_Survdata$pack_years_calculated +  
                       Tx321_tumour_Survdata$adjuvant_treatment_YN +
                       Tx321_tumour_Survdata$TNM_combine +
                       Tx321_tumour_Survdata$IASLC_LUAD_grade
)

MVA.cox.res <- summary(MVA.cox.res)
MVA.cox.res <- data.frame(HR = MVA.cox.res$coefficients[,2], lower_ci = MVA.cox.res$conf.int[,3], upper_ci = MVA.cox.res$conf.int[,4], P = MVA.cox.res$coefficients[,5])
MVA.cox.res$Predictors <- c("Mean RS","Female","Age","pack_years","Adjuvant treatment","Stage II","Stage III","Grade 2","Grade 3")

#add ref
idx <- min(grep("Female", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","Male")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Adjuvant", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","No adjuvant")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Stage", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","Stage I")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Grade", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","IMA")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

#turn data from chr to numeric
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#tidy up presenting values - round up and 95% CI
MVA.cox.res$pvalue <-signif(MVA.cox.res$P,digits = 1) 

#add haeder
MVA.cox.res <- rbind(c(1,NA,NA,"","","p-value","HR","CI"),MVA.cox.res)
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#Forest plot
forestplot(MVA.cox.res[,c(5:6)],graph.pos=2,
           mean = MVA.cox.res$HR,
           lower = MVA.cox.res$lower_ci,
           upper = MVA.cox.res$upper_ci,
           xlog=TRUE,title="Hazard Ratio",
           boxsize = 0.4,ci.vertices=TRUE,ci.vertices.height = 0.2,
           txt_gp=fpTxtGp(label=gpar(cex=0.9),
                          ticks=gpar(cex=0.9),
                          xlab=gpar(cex = 1),
                          title=gpar(cex = 1.2)),lwd.ci=1,colgap=unit(6,"mm"),
           col=fpColors(box="black", lines="black", zero = "gray50"),
           graphwidth = unit(12,"cm"))


## 2C

sample.df  <- Tx321_riskscore
colnames(sample.df)[5] <- "ORACLE_RS"

#iterate single sampling 1000 times
sigs <- c("ORACLE","Li","Wang","Zhao","Song","Jin","LiFeng")

#ORACLE
ORACLE.surv.df <- list()
set.seed(132)
for (i in 1:1000){
  ORACLE.surv.df[[i]] <- create_surv(sigs[1])
}  

#Li
Li.surv.df <- list()
set.seed(133)
for (i in 1:1000){
  Li.surv.df[[i]] <- create_surv(sigs[2])
}  

#Wang
Wang.surv.df <- list()
set.seed(134)
for (i in 1:1000){
  Wang.surv.df[[i]] <- create_surv(sigs[3])
}  

#Zhao
Zhao.surv.df <- list()
set.seed(135)
for (i in 1:1000){
  Zhao.surv.df[[i]] <- create_surv(sigs[4])
}  

#Song
Song.surv.df <- list()
set.seed(136)
for (i in 1:1000){
  Song.surv.df[[i]] <- create_surv(sigs[5])
}  

#Jin
Jin.surv.df <- list()
set.seed(137)
for (i in 1:1000){
  Jin.surv.df[[i]] <- create_surv(sigs[6])
}  

#LiFeng
LiFeng.surv.df <- list()
set.seed(138)
for (i in 1:1000){
  LiFeng.surv.df[[i]] <- create_surv(sigs[7])
}  

#data frame for cox regression results
#ORACLE
ORACLE.cox.df <- list()
for(i in 1:1000){
  ORACLE.cox.df[[i]] <- bootsrap_surv(survdf = ORACLE.surv.df[[i]])
}

#Li
Li.cox.df <- list()
for(i in 1:1000){
  Li.cox.df[[i]] <- bootsrap_surv(survdf = Li.surv.df[[i]])
}

#Wang
Wang.cox.df <- list()
for(i in 1:1000){
  Wang.cox.df[[i]] <- bootsrap_surv(survdf = Wang.surv.df[[i]])
}

#Zhao
Zhao.cox.df <- list()
for(i in 1:1000){
  Zhao.cox.df[[i]] <- bootsrap_surv(survdf = Zhao.surv.df[[i]])
}

#Song
Song.cox.df <- list()
for(i in 1:1000){
  Song.cox.df[[i]] <- bootsrap_surv(survdf = Song.surv.df[[i]])
}

#Jin
Jin.cox.df <- list()
for(i in 1:1000){
  Jin.cox.df[[i]] <- bootsrap_surv(survdf = Jin.surv.df[[i]])
}

#LiFeng
LiFeng.cox.df <- list()
for(i in 1:1000){
  LiFeng.cox.df[[i]] <- bootsrap_surv(survdf = LiFeng.surv.df[[i]])
}

#tidy up - list to data frame
ORACLE.cox.df <- Reduce(rbind,ORACLE.cox.df)
Li.cox.df <- Reduce(rbind,Li.cox.df)
Wang.cox.df <- Reduce(rbind,Wang.cox.df)
Zhao.cox.df <- Reduce(rbind,Zhao.cox.df)
Song.cox.df <- Reduce(rbind,Song.cox.df)
Jin.cox.df <- Reduce(rbind,Jin.cox.df)
LiFeng.cox.df <- Reduce(rbind,LiFeng.cox.df)

#calculate boot HR
ORACLE_HR_mu <- mean(ORACLE.cox.df$HR)
ORACLE_HR_boot_lci <- mean(ORACLE.cox.df$cox_lci) #lower bound
ORACLE_HR_boot_uci <- mean(ORACLE.cox.df$cox_uci) #higher bound

#
Li_HR_mu <- mean(Li.cox.df$HR)
Li_HR_boot_lci <- mean(Li.cox.df$cox_lci) #lower bound
Li_HR_boot_uci <- mean(Li.cox.df$cox_uci) #higher bound

#
Wang_HR_mu <- mean(Wang.cox.df$HR)
Wang_HR_boot_lci <- mean(Wang.cox.df$cox_lci) #lower bound
Wang_HR_boot_uci <- mean(Wang.cox.df$cox_uci) #higher bound

#
Zhao_HR_mu <- mean(Zhao.cox.df$HR)
Zhao_HR_boot_lci <- mean(Zhao.cox.df$cox_lci) #lower bound
Zhao_HR_boot_uci <- mean(Zhao.cox.df$cox_uci) #higher bound

#
Song_HR_mu <- mean(Song.cox.df$HR)
Song_HR_boot_lci <- mean(Song.cox.df$cox_lci) #lower bound
Song_HR_boot_uci <- mean(Song.cox.df$cox_uci) #higher bound

#
Jin_HR_mu <- mean(Jin.cox.df$HR)
Jin_HR_boot_lci <- mean(Jin.cox.df$cox_lci) #lower bound
Jin_HR_boot_uci <- mean(Jin.cox.df$cox_uci) #higher bound


#
LiFeng_HR_mu <- mean(LiFeng.cox.df$HR)
LiFeng_HR_boot_lci <- mean(LiFeng.cox.df$cox_lci) #lower bound
LiFeng_HR_boot_uci <- mean(LiFeng.cox.df$cox_uci) #higher bound

#combine data frame
ORACLE.cox.df$cohort <- "ORACLE"
Li.cox.df$cohort <- "Li"
Wang.cox.df$cohort <- "Wang"
Zhao.cox.df$cohort <- "Zhao"
Song.cox.df$cohort <- "Song"
Jin.cox.df$cohort <- "Jin"
LiFeng.cox.df$cohort <- "LiFeng"

ggplot(data=ORACLE.cox.df,aes(HR))+
  geom_ribbon(aes(xmax = mean(cox_lci), xmin = mean(cox_uci),ymax = 3, ymin= 0),alpha=0.2, fill = "#FDAE5F") + 
  geom_density(fill="#4EB0EC",alpha=0.4)+ 
  xlab("Hazard ratio") + ylab("Density")+
  scale_x_continuous(breaks = seq(0,3.5,1),expand = c(0,0),limits = c(0,3.5))+
  scale_y_continuous(breaks = seq(0,3,1),expand = c(0,0),limits = c(0,3))+
  geom_vline(xintercept = 1, lty="dashed",col="#25508F")+
  geom_vline(xintercept = ORACLE_HR_mu, col = "red") +
  theme(plot.margin = unit(c(2,4,2,2),"mm"),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))


## 2D

Survdata_stageI <- Tx321_tumour_Survdata[which(Tx321_tumour_Survdata$pTNMStage_v8 == "1a" | Tx321_tumour_Survdata$pTNMStage_v8 == "1b"),]

#subgroup - substaging
j <- ggsurvplot(fit = survfit(Surv(os_time,cens_os)~pTNMStage_v8,data = Survdata_stageI),
                legend.labs = c("IA", "IB"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","#EE0000FF"),size = 1.5, 
                title="Substaging - StageI", 
                legend="none", xlab="Time (Years)", ylab="Overall Survival",
                pval = T,pval.size = 8,
                ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j
coxmd <- coxph(Surv(os_time, cens_os) ~ pTNMStage_v8 , data = Survdata_stageI)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])

#subgroup - ORACLE
Survdata_stageI$ORACLE_class <- factor(Survdata_stageI$ORACLE_class,levels = c("Low","Discordant","High"))
j<-ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Survdata_stageI),
              legend.labs = c("Low", "Discordant","High"),
              risk.table = TRUE, 
              palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
              title="ORACLE", 
              legend="none", xlab="Time (Years)", ylab="Overall Survival",
              pval = T,pval.size = 8,
              ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                              panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
              tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j
coxmd <- coxph(Surv(os_time, cens_os) ~ ORACLE_class , data = Survdata_stageI)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])




### Figure 3

## 3A

j <- ggsurvplot(fit = survfit(Surv(lung_specific_time,cens_lung_specific)~ORACLE_class,data = Tx321_tumour_Survdata),
                legend.labs = c( "Low","Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                legend="none", xlab="Time (years)", ylab="Lung-specific Survival",
                pval = T,pval.size = 6,pval.coord = c(0, 0.05),
                ggtheme = theme(title = element_text(size = 18),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j


## 3B

j <- ggsurvplot(fit = survfit(Surv(lung_specific_time,cens_lung_specific)~ORACLE_class,data = Survdata_stageI),
                legend.labs = c( "Low","Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                legend="none", xlab="Time (years)", ylab="Lung-specific Survival",
                pval = T,pval.size = 6,pval.coord = c(0, 0.05),
                ggtheme = theme(title = element_text(size = 18),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j


## 3C

#calculate ORACLE risk-score
Mascaux_ORACLE_RS <- data.frame(riskscore = colSums(Mascaux_ORACLE[Supp_Table_5_ORACLE$Gene.Symbol,]*Supp_Table_5_ORACLE$Model.Coefficient))
rownames(Mascaux_ORACLE_RS) <- colnames(patient.annot)[-1]
Mascaux_ORACLE_RS$patient <- gsub("-.*", rownames(Mascaux_ORACLE_RS),replacement = "")

#development group
Mascaux_ORACLE_RS$group <- NA
Mascaux_ORACLE_RS$group[which(grepl("hyperplasia",rownames(Mascaux_ORACLE_RS)))] <- "hyperplasia"
Mascaux_ORACLE_RS$group[which(grepl("normal",rownames(Mascaux_ORACLE_RS)))] <- "normal"
Mascaux_ORACLE_RS$group[which(grepl("metaplasia",rownames(Mascaux_ORACLE_RS)))] <- "metaplasia"
Mascaux_ORACLE_RS$group[which(grepl("mild",rownames(Mascaux_ORACLE_RS)))] <- "mild\ndysplasia"
Mascaux_ORACLE_RS$group[which(grepl("moderate",rownames(Mascaux_ORACLE_RS)))] <- "moderate\ndysplasia"
Mascaux_ORACLE_RS$group[which(grepl("severe",rownames(Mascaux_ORACLE_RS)))] <- "severe\ndysplasia"
Mascaux_ORACLE_RS$group[which(grepl("carcinoma.in.situ",rownames(Mascaux_ORACLE_RS)))] <- "CIS"
Mascaux_ORACLE_RS$group[which(grepl("squamous",rownames(Mascaux_ORACLE_RS)))] <- "SCC"

#grade by Mascaux
Mascaux_ORACLE_RS$grade <- NA
Mascaux_ORACLE_RS$grade[which(grepl("normofluorescent|hypofluorescent|hyperplasia",rownames(Mascaux_ORACLE_RS)))] <- "Normal"
Mascaux_ORACLE_RS$grade[which(grepl("metaplasia|mild|moderate",rownames(Mascaux_ORACLE_RS)))] <- "Low-grade"
Mascaux_ORACLE_RS$grade[which(grepl("severe|carcinoma.in.situ",rownames(Mascaux_ORACLE_RS)))] <- "High-grade"
Mascaux_ORACLE_RS$grade[which(grepl("squamous",rownames(Mascaux_ORACLE_RS)))] <- "SCC"

ggplot(Mascaux_ORACLE_RS,aes(group,riskscore))+geom_boxplot(aes(fill=group))+
  xlab("Developmental stages")+ylab("ORACLE Risk score")+
  scale_y_continuous(expand = c(0,0),limits = c(0,3))+
  scale_fill_manual(values = brewer.pal(8,"Oranges"))+
  geom_text(data = r ,aes(x = group, y = 2.8, label = signif(pval,digits = 2))) +
  theme(legend.position="none",panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))


## 3D

#paired analysis on whole LUAD cohort
Rec.sample.df <- Tx421_riskscore
Rec.sample.df$PrimaryRegion <- gsub("-", Rec.sample.df$sample_name_cruk, replacement = ".")

#seeding/non-seeding info
Rec.sample.df <- left_join(Rec.sample.df, seed.region[,c("PrimaryRegion","Metastasizing")], by ="PrimaryRegion")
Rec.sample.df <- Rec.sample.df[which(!is.na(Rec.sample.df$Metastasizing)),]

#keep samples with paired seeding and one non-seeding
remove <- table(Rec.sample.df$cruk_id, Rec.sample.df$Metastasizing) %>% as.data.frame()
Rec.sample.df <- Rec.sample.df[-which(Rec.sample.df$cruk_id %in% remove$Var1[which(remove$Freq == 0)]),]

#select recurrence samples
Rec.sample.df <- left_join(Rec.sample.df, TRACERx_clinical[,c("cruk_id","Relapse_cat_new","Relapse_cat")], by = "cruk_id")
Rec.sample.df <- Rec.sample.df[which(Rec.sample.df$Relapse_cat_new != "No rec"),]

#assign met-seed/non-met-seed
Rec.sample.df$Metastasizing[which(Rec.sample.df$Metastasizing == "TRUE")] <- "Seeding"
Rec.sample.df$Metastasizing[which(Rec.sample.df$Metastasizing == "FALSE")] <- "Non-seeding"

#linear mixed effect model for multiple samples in a patient
md <-nlme::lme(data=Rec.sample.df,ORACLE_riskscore~Metastasizing,na.action = na.omit,random=~1|cruk_id)
r <- summary(md)

ggplot(Rec.sample.df,aes(Metastasizing,ORACLE_riskscore))+
  geom_boxplot(aes(fill=Metastasizing), width=0.5)+
  geom_jitter(width = 0.2)+
  ylab("ORACLE Risk-score")+xlab("")+
  scale_fill_manual(values=c("#7DB5FF","#FDCF0B","#FF899F"))+
  scale_y_continuous(breaks=seq(9,12,1),expand = c(0,0),limits = c(9,12))+
  geom_text(aes(x=2,y=11.8,label=paste0("p = ",signif(r$tTable[2,5],digits = 2))),size=4) +
  theme(legend.position = "none",panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())


## 3E

Relapse_survdata$ORACLE_class <- factor(Relapse_survdata$ORACLE_class, levels = c("Low","Discordant","High"))
j <- ggsurvplot(fit = survfit(Surv(dfs_time,cens_dfs)~ORACLE_class,data = Relapse_survdata),
                legend.labs = c( "Low","Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                legend="none", xlab="Time (years)", ylab="Disease-Free Survival",
                pval = T,pval.size = 6,pval.coord = c(0, 0.05),
                ggtheme = theme(title = element_text(size = 18),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j
coxmd <- coxph(Surv(dfs_time, cens_dfs) ~ ORACLE_class , data = Relapse_survdata)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])


## 3F

Relapse_survdata_stageI <- Relapse_survdata[which(Relapse_survdata$pTNMStage_v8 %in% c("1a","1b")),]
j <- ggsurvplot(fit = survfit(Surv(dfs_time,cens_dfs)~ORACLE_class,data = Relapse_survdata_stageI),
                legend.labs = c( "Low","Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                legend="none", xlab="Time (years)", ylab="Disease-Free Survival",
                pval = T,pval.size = 6,pval.coord = c(0, 0.05),
                ggtheme = theme(title = element_text(size = 18),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j
coxmd <- coxph(Surv(dfs_time, cens_dfs) ~ ORACLE_class , data = Relapse_survdata)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])




### Figure 4

## 4A
IC50_results$cor <- ifelse(IC50_results$estimate > 0 , "Positive","Negative")

#label FDA drugs for NSCLC
IC50_results$label <- IC50_results$cor
IC50_results$label[which(IC50_results$significance=="FALSE")] <- NA
IC50_results$label[which(IC50_results$FDA=="TRUE")] <- "FDA"

ggplot(IC50_results,aes(x = estimate,y = -log10(pvalue)))+
  geom_point(alpha=0.7,aes(col=significance),size=3)+
  geom_point(data = IC50_results[which(IC50_results$FDA=="TRUE"),],alpha=0.8,size=3,aes(x = estimate,y = -log10(pvalue)),pch=21,col="black",fill="#D1D1D1")+
  geom_point(data = IC50_results[which(IC50_results$Name=="cisplatin"),],alpha=0.8,size=3,aes(x = estimate,y = -log10(pvalue)),pch=21,col="black",fill="#E20000")+
  scale_x_continuous(breaks=seq(-0.5,0.5,0.5),expand=c(0,0),limits = c(-0.5,0.5))+
  ylab("-log(pvalue)") + xlab("") +
  ggtitle("")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  geom_text_repel(data = IC50_results[which(IC50_results$FDA=="TRUE"),],aes(label = Name),min.segment.length = unit(0.1, "lines"), size=3,max.overlaps = 20)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = -log10(0.05),lty="dashed") + 
  scale_color_manual(values=c("#D1D1D1","#E20000"),na.value = "#D1D1D1")+
  theme(legend.text = element_text(size=12),axis.title = element_text(size=14),axis.text = element_text(size = 12),panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())


## 4B

ggplot(IC50_results, aes(fct_reorder(pathway,estimate),estimate)) +  
  geom_point(aes(col=significance),size=3,alpha=0.8) + 
  geom_boxplot(width=0.4,col="#0060A2",alpha=0,outlier.colour = NULL) +
  scale_color_manual(values = c("azure4","#E20000")) + 
  scale_y_continuous(expand = c(0,0),limits = c(-0.5,0.5)) +
  ylab("") + xlab("Targeting pathway") + 
  geom_hline(yintercept = 0, lty="dashed") +
  theme_bw() +
  theme(legend.position = "none",axis.text.x = element_text(angle=90, vjust=0.5,hjust = 1),panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())


## 4C

#split into adjuvant/no adjuvant cohorts
Survdata_No_adjuvant <- Tx321_tumour_Survdata[which( Tx321_tumour_Survdata$adjuvant_treatment_YN == "No adjuvant"),]
Survdata_with_adjuvant <- Tx321_tumour_Survdata[which( Tx321_tumour_Survdata$adjuvant_treatment_YN == "Adjuvant"),]

#No adjuvant
j <- ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Survdata_No_adjuvant),
                legend.labs = c("Low", "Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                title="No adjuvant", 
                legend="none", xlab="Time (Years)", ylab="Overall Survival",
                pval = T,pval.size = 8,
                ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j
coxmd <- coxph(Surv(os_time, cens_os) ~ ORACLE_class , data = Survdata_No_adjuvant)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])

#With adjuvant
j<-ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Survdata_with_adjuvant),
              legend.labs = c("Low", "Discordant","High"),
              risk.table = TRUE, 
              palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
              title="With adjuvant", 
              legend="none", xlab="Time (Years)", ylab="Overall Survival",
              pval = T,pval.size = 8,
              ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                              panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
              tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j
coxmd <- coxph(Surv(os_time, cens_os) ~ ORACLE_class , data = Survdata_with_adjuvant)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])





### Figure 5

## 5A

#MLR with clinical model
MLR_model_clinical <- lm(data=Clinical.scale,mean~age+sexFemale+`smoking_statusEx-Smoker`+`smoking_statusCurrent Smoker`+ `pTNMStage_v8II`+ `pTNMStage_v8III` + biopsy_no+Volume+KI67_score)
r <- summary(MLR_model_clinical)
pvalue.MLR.clinical <- data.frame(correlate = names(r$coefficients[-1,1]),pvalue = r$coefficients[-1,4],coef = r$coefficients[-1,1])
pvalue.MLR.clinical$correlate <- c("Age","Female","Ex-Smoker","Smoker","Stage II","Stage III","#biopsy","Volume","Ki67")
pvalue.MLR.clinical$correlate <- factor(pvalue.MLR.clinical$correlate , levels = pvalue.MLR.clinical$correlate)

pvalue.MLR.clinical$significance <- ""
pvalue.MLR.clinical$significance[which(pvalue.MLR.clinical$pvalue < 0.05 & pvalue.MLR.clinical$pvalue >= 0.01)] <- "*"
pvalue.MLR.clinical$significance[which(pvalue.MLR.clinical$pvalue < 0.01 & pvalue.MLR.clinical$pvalue >= 0.005)] <- "**"
pvalue.MLR.clinical$significance[which(pvalue.MLR.clinical$pvalue < 0.005)] <- "***"

#scale genetic metric values
Genetics.scale <-  dplyr::select(Features.df, cruk_id, mean, frac_abberant_genom_subcl, wFLOH, num_clonal_gds, num_subclonal_gds, num_clonal_drivers, num_subclonal_drivers, RecentSubclExpansionScore)
Genetics.scale[,-1] <- scale(Genetics.scale[,-1], center = T)

#MLR with genetics model
MLR_model_genetic <- lm(data = Genetics.scale, mean ~  frac_abberant_genom_subcl + wFLOH + num_clonal_gds + num_subclonal_gds + num_clonal_drivers + num_subclonal_drivers + RecentSubclExpansionScore  ) 
r <- summary(MLR_model_genetic)
pvalue.MLR.genetics <- data.frame(correlate = names(r$coefficients[-1,1]),pvalue = r$coefficients[-1,4],coef = r$coefficients[-1,1])
pvalue.MLR.genetics$correlate <- c("SCNA-ITH","FLOH","Clonal WGD","Subclonal WGD","Clonal drivers","Subclonal drivers","Recent subclonal expansion")
pvalue.MLR.genetics$correlate <- factor(pvalue.MLR.genetics$correlate , levels = pvalue.MLR.genetics$correlate)

pvalue.MLR.genetics$significance <- ""
pvalue.MLR.genetics$significance[which(pvalue.MLR.genetics$pvalue < 0.1 & pvalue.MLR.genetics$pvalue >= 0.01)] <- "*"
pvalue.MLR.genetics$significance[which(pvalue.MLR.genetics$pvalue < 0.01 & pvalue.MLR.genetics$pvalue >= 0.005)] <- "**"
pvalue.MLR.genetics$significance[which(pvalue.MLR.genetics$pvalue < 0.005)] <- "***"

MLR.results <- rbind(pvalue.MLR.clinical,pvalue.MLR.genetics)

ggplot(MLR.results)+
  geom_tile(aes(correlate, y = 1,fill=coef))+
  scale_y_continuous(expand = c(0,0))+
  xlab("")+ylab("Mean")+labs(col="")+
  scale_fill_distiller(name = '',limits=c(-0.45,0.45),palette ="RdBu",na.value = "white") +
  scale_x_discrete(expand = c(0,0)) +
  geom_text(aes(x = correlate, y = 1, label=significance),size=5,hjust=0.5)+
  theme(legend.position = "top",axis.text = element_blank(), axis.title.y = element_text(angle = 0, vjust=0.5,hjust = 1),axis.ticks = element_blank(),panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())


## 5B

#level
Tx_biomarkers$adjuvant_treatment_YN <- factor(Tx_biomarkers$adjuvant_treatment_YN, levels = c("No adjuvant","Adjuvant"))
Tx_biomarkers$pTNMStage_v8 <- factor(Tx_biomarkers$pTNMStage_v8, levels = c("I","II","III"))
Tx_biomarkers$sex <- factor(Tx_biomarkers$sex, levels = c("Male","Female"))

MVA.cox.res <- coxph(data = Tx_biomarkers,formula = Surv(Tx_biomarkers$os_time,Tx_biomarkers$cens_os) ~ 
                       age +
                       sex +
                       pTNMStage_v8 +
                       pack_years_calculated +
                       adjuvant_treatment_YN +
                       SCNA_ITH + 
                       as.numeric(subclonal_wgd) +
                       RecentSubclExpansionScore +
                       as.numeric(ctDNA_status) +
                       as.numeric(STAS) + 
                       ORACLE_mean)

ggforest(data = as.data.frame(Tx_biomarkers), MVA.cox.res)


## 5C

#scale ORACLE score
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
Tx_biomarker_avail$ORACLE_scaled <- range01(Tx_biomarker_avail$ORACLE_mean)

#iterate from FU time = 1 to 5
glm_time_all <- list()
for(i in 1:5){
  glm_time_all[[i]] <- biomarker.glm(i)
}
glm_time_all <- Reduce(rbind,glm_time_all)

#ordered by the first year PVE
glm_time_all <- dplyr::arrange(glm_time_all, OS_time, -Explained_var)
glm_time_all$variable <- factor(glm_time_all$variable, levels = glm_time_all$variable[length(unique(glm_time_all$variable)):1])

#bar plot
ggplot(glm_time_all, aes(OS_time, Explained_var)) + xlab("Time (years)") + ylab("% variance explained") + 
  geom_col(aes(fill=variable),col="black",alpha=0.8) +
  scale_fill_manual(values = brewer.pal(12,"Paired")[c(6,4,9,3,7,1)]) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,40)) +
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))





### Extended Data Fig 2

## ED 2A
Li_plot <- ggplot(signature_TSB,aes(x=fct_reorder(cruk_id, as.numeric(Li_class) + Li_RS, .fun = mean),y = Li_RS))+
  geom_point(aes(col = Li_class),alpha = 0.5,pch=16,size = 3)+
  geom_hline(yintercept = median(signature_TSB$Li_RS),lty = "dotted")+
  geom_line(col = "black")+
  xlab("")+
  ylab("Li Risk score")+
  labs(color="Class")+
  scale_y_continuous(breaks = seq(-0.5,0.5,0.2),expand = c(0,0),limits = c(-0.55,0.55))+
  scale_color_manual(values = c("#3B4992FF","azure4","#EE0000FF"))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1,size=6),axis.ticks.x = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

Wang_plot <- ggplot(signature_TSB,aes(x=fct_reorder(cruk_id, as.numeric(Wang_class) + Wang_RS, .fun = mean),y = Wang_RS))+
  geom_point(aes(col = Wang_class),alpha = 0.5,pch=16,size = 3)+
  geom_hline(yintercept = median(signature_TSB$Wang_RS),lty = "dotted")+
  geom_line(col = "black")+
  xlab("")+
  ylab("Wang Risk score")+
  labs(color="Class")+
  scale_y_continuous(breaks = seq(0,2,0.4),expand = c(0,0),limits = c(0,2))+
  scale_color_manual(values = c("#3B4992FF","azure4","#EE0000FF"))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1,size=6),axis.ticks.x = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

Zhao_plot <- ggplot(signature_TSB,aes(x=fct_reorder(cruk_id, as.numeric(Zhao_class) + Zhao_RS, .fun = mean),y = Zhao_RS))+
  geom_point(aes(col = Zhao_class),alpha = 0.5,pch=16,size = 3)+
  geom_hline(yintercept = median(signature_TSB$Zhao_RS),lty = "dotted")+
  geom_line(col = "black")+
  xlab("")+
  ylab("Zhao Risk score")+
  labs(color="Class")+
  scale_y_continuous(breaks = seq(0,75,25),expand = c(0,0),limits = c(-4.5,75))+
  scale_color_manual(values = c("#3B4992FF","azure4","#EE0000FF"))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1,size=6),axis.ticks.x = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

Song_plot <- ggplot(signature_TSB,aes(x=fct_reorder(cruk_id, as.numeric(Song_class) + Song_RS, .fun = mean),y = Song_RS))+
  geom_point(aes(col = Song_class),alpha = 0.5,pch=16,size = 3)+
  geom_hline(yintercept = median(signature_TSB$Song_RS),lty = "dotted")+
  geom_line(col = "black")+
  xlab("")+
  ylab("Song Risk score")+
  labs(color="Class")+
  scale_y_continuous(breaks = seq(0,6,2),expand = c(0,0),limits = c(0,6))+
  scale_color_manual(values = c("#3B4992FF","azure4","#EE0000FF"))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1,size=6),axis.ticks.x = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

Jin_plot <- ggplot(signature_TSB,aes(x=fct_reorder(cruk_id, as.numeric(Jin_class) + Jin_RS, .fun = mean),y = Jin_RS))+
  geom_point(aes(col = Jin_class),alpha = 0.5,pch=16,size = 3)+
  geom_hline(yintercept = median(signature_TSB$Jin_RS),lty = "dotted")+
  geom_line(col = "black")+
  xlab("")+
  ylab("Jin Risk score")+
  labs(color="Class")+
  scale_y_continuous(breaks = seq(-8,-3,1),expand = c(0,0),limits = c(-8,-3))+
  scale_color_manual(values = c("#3B4992FF","azure4","#EE0000FF"))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1,size=6),axis.ticks.x = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

LiFeng_plot <- ggplot(signature_TSB,aes(x=fct_reorder(cruk_id, as.numeric(LiFeng_class) + LiFeng_RS, .fun = mean),y = LiFeng_RS))+
  geom_point(aes(col = LiFeng_class),alpha = 0.5,pch=16,size = 3)+
  geom_hline(yintercept = median(signature_TSB$LiFeng_RS),lty = "dotted")+
  geom_line(col = "black")+
  xlab("")+
  ylab("LiFeng Risk score")+
  labs(color="Class")+
  scale_y_continuous(breaks = seq(-3,3,2),expand = c(0,0),limits = c(-3,3))+
  scale_color_manual(values = c("#3B4992FF","azure4","#EE0000FF"))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1,size=6),axis.ticks.x = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#grob
Li_plot <- ggplotGrob(Li_plot + theme(plot.margin = unit(c(0, 4, 0, 2), "mm")))
Wang_plot <- ggplotGrob(Wang_plot + theme(plot.margin = unit(c(0, 4, 0, 2), "mm")))
Zhao_plot <- ggplotGrob(Zhao_plot + theme(plot.margin = unit(c(0, 4, 0, 2), "mm")))
Song_plot <- ggplotGrob(Song_plot + theme(plot.margin = unit(c(0, 4, 0, 2), "mm")))
Jin_plot <- ggplotGrob(Jin_plot + theme(plot.margin = unit(c(0, 4, 0, 2), "mm")))
LiFeng_plot <- ggplotGrob(LiFeng_plot + theme(plot.margin = unit(c(0, 4, 2, 2), "mm")))

#adjust height and combine plots
g <- gtable_rbind(Li_plot,Wang_plot,Zhao_plot,Song_plot,Jin_plot,LiFeng_plot)
id_panels_h <- unique(g$layout[g$layout$name=="panel","t"])
g$heights[id_panels_h] <- grid::unit(c(4,4,4,4,4,4),"null")
grid.arrange(g)


## ED 2B

#join stage info
TSB_summary <- left_join(TSB_summary,all_tumour_df_20220209[,c("tumour_id_mphase_cruk","pTNMStage_v8")], by = "tumour_id_mphase_cruk")
TSB_summary$pTNMStage_v8[which(TSB_summary$cruk_id == "CRUK0704")] <- "2b" #pathology report - T3N0M0

#classify into stage I, II, III
TSB_summary$pTNMStage_v8[which(TSB_summary$pTNMStage_v8 %in% c("1a","1b"))] <- "Stage I"
TSB_summary$pTNMStage_v8[which(TSB_summary$pTNMStage_v8 %in% c("2a","2b"))] <- "Stage II"
TSB_summary$pTNMStage_v8[which(TSB_summary$pTNMStage_v8 %in% c("3a","3b"))] <- "Stage III"

#level risk-class
TSB_summary <- dplyr::select(TSB_summary, cruk_id, ORACLE_class, Li_class, Wang_class, Zhao_class, Song_class, Jin_class, LiFeng_class, pTNMStage_v8)
TSB_summary$ORACLE_class <- factor(TSB_summary$ORACLE_class, levels = c("Low","High","Discordant"))
TSB_summary$Li_class <- factor(TSB_summary$Li_class, levels = c("Low","High","Discordant"))
TSB_summary$Wang_class <- factor(TSB_summary$Wang_class, levels = c("Low","High","Discordant"))
TSB_summary$Zhao_class <- factor(TSB_summary$Zhao_class, levels = c("Low","High","Discordant"))
TSB_summary$Song_class <- factor(TSB_summary$Song_class, levels = c("Low","High","Discordant"))
TSB_summary$Jin_class <- factor(TSB_summary$Jin_class, levels = c("Low","High","Discordant"))
TSB_summary$LiFeng_class <- factor(TSB_summary$LiFeng_class, levels = c("Low","High","Discordant"))

#plot risk-class ~ stage for each signature

#ORACLE risk-class frequency
df <- expand.grid( Risk = c("Low","High","Discordant"),Stage = c("Stage I","Stage II","Stage III"), perc = NA)
df$perc <- rep(table(TSB_summary$pTNMStage_v8)*100/nrow(TSB_summary),each=3)
df$risk.perc <- c(table(TSB_summary$ORACLE_class,TSB_summary$pTNMStage_v8)[,1]/sum(table(TSB_summary$ORACLE_class,TSB_summary$pTNMStage_v8)[,1]),
                  table(TSB_summary$ORACLE_class,TSB_summary$pTNMStage_v8)[,2]/sum(table(TSB_summary$ORACLE_class,TSB_summary$pTNMStage_v8)[,2]),
                  table(TSB_summary$ORACLE_class,TSB_summary$pTNMStage_v8)[,3]/sum(table(TSB_summary$ORACLE_class,TSB_summary$pTNMStage_v8)[,3]))

#chi-squared test
chi.table <- TSB_summary
chi.table$ORACLE_class <- as.character(chi.table$ORACLE_class)
chi.table$ORACLE_class[which(chi.table$ORACLE_class %in% c("High","Low"))] <- "Non-discordant"
chisq_p <- chisq.test(table(chi.table$ORACLE_class,chi.table$pTNMStage_v8))$p.value
ggplot(df,aes(Stage,risk.perc*100,fill=Risk)) + geom_col(width=0.8,position = "dodge") +
  ylab("Survival risk classification (%)")+
  xlab("")+ggtitle("ORACLE")+
  geom_text(aes(y = risk.perc*100+1,label = paste0(signif(risk.perc*100,digits = 2),"%")),position = position_dodge(width=0.8),size=3)+
  geom_text(data=df, aes(x=2.5,y=75,label = paste0("chi-sq test, p=",signif(chisq_p,digits = 1))),size=5)+
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0, 0), limits = c(0,100))+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme(legend.position = "none",axis.text.x = element_text(size=12),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#Li risk-class frequency
df <- expand.grid( Risk = c("Low","High","Discordant"),Stage = c("Stage I","Stage II","Stage III"), perc = NA)
df$perc <- rep(table(TSB_summary$pTNMStage_v8)*100/nrow(TSB_summary),each=3)
df$risk.perc <- c(table(TSB_summary$Li_class,TSB_summary$pTNMStage_v8)[,1]/sum(table(TSB_summary$Li_class,TSB_summary$pTNMStage_v8)[,1]),
                  table(TSB_summary$Li_class,TSB_summary$pTNMStage_v8)[,2]/sum(table(TSB_summary$Li_class,TSB_summary$pTNMStage_v8)[,2]),
                  table(TSB_summary$Li_class,TSB_summary$pTNMStage_v8)[,3]/sum(table(TSB_summary$Li_class,TSB_summary$pTNMStage_v8)[,3]))

#chi-squared test
chi.table <- TSB_summary
chi.table$Li_class <- as.character(chi.table$Li_class)
chi.table$Li_class[which(chi.table$Li_class %in% c("High","Low"))] <- "Non-discordant"
chisq_p <- chisq.test(table(chi.table$Li_class,chi.table$pTNMStage_v8))$p.value

ggplot(df,aes(Stage,risk.perc*100,fill=Risk)) + geom_col(width=0.8,position = "dodge") +
  ylab("Survival risk classification (%)")+
  xlab("")+ggtitle("Li")+
  geom_text(aes(y = risk.perc*100+1,label = paste0(signif(risk.perc*100,digits = 2),"%")),position = position_dodge(width=0.8),size=3)+
  geom_text(data=df, aes(x=2.5,y=75,label = paste0("chi-sq test, p=",signif(chisq_p,digits = 1))),size=5)+
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0, 0), limits = c(0,100))+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme(legend.position = "none",axis.text.x = element_text(size=12),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#Wang risk-class frequency
df <- expand.grid( Risk = c("Low","High","Discordant"),Stage = c("Stage I","Stage II","Stage III"), perc = NA)
df$perc <- rep(table(TSB_summary$pTNMStage_v8)*100/nrow(TSB_summary),each=3)
df$risk.perc <- c(table(TSB_summary$Wang_class,TSB_summary$pTNMStage_v8)[,1]/sum(table(TSB_summary$Wang_class,TSB_summary$pTNMStage_v8)[,1]),
                  table(TSB_summary$Wang_class,TSB_summary$pTNMStage_v8)[,2]/sum(table(TSB_summary$Wang_class,TSB_summary$pTNMStage_v8)[,2]),
                  table(TSB_summary$Wang_class,TSB_summary$pTNMStage_v8)[,3]/sum(table(TSB_summary$Wang_class,TSB_summary$pTNMStage_v8)[,3]))

#chi-squared test
chi.table <- TSB_summary
chi.table$Wang_class <- as.character(chi.table$Wang_class)
chi.table$Wang_class[which(chi.table$Wang_class %in% c("High","Low"))] <- "Non-discordant"
chisq_p <- chisq.test(table(chi.table$Wang_class,chi.table$pTNMStage_v8))$p.value

ggplot(df,aes(Stage,risk.perc*100,fill=Risk)) + geom_col(width=0.8,position = "dodge") +
  ylab("Survival risk classification (%)")+
  xlab("")+ggtitle("Wang")+
  geom_text(aes(y = risk.perc*100+1,label = paste0(signif(risk.perc*100,digits = 2),"%")),position = position_dodge(width=0.8),size=3)+
  geom_text(data=df, aes(x=2.5,y=75,label = paste0("chi-sq test, p=",signif(chisq_p,digits = 1))),size=5)+
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0, 0), limits = c(0,100))+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme(legend.position = "none",axis.text.x = element_text(size=12),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#Zhao risk-class frequency
df <- expand.grid( Risk = c("Low","High","Discordant"),Stage = c("Stage I","Stage II","Stage III"), perc = NA)
df$perc <- rep(table(TSB_summary$pTNMStage_v8)*100/nrow(TSB_summary),each=3)
df$risk.perc <- c(table(TSB_summary$Zhao_class,TSB_summary$pTNMStage_v8)[,1]/sum(table(TSB_summary$Zhao_class,TSB_summary$pTNMStage_v8)[,1]),
                  table(TSB_summary$Zhao_class,TSB_summary$pTNMStage_v8)[,2]/sum(table(TSB_summary$Zhao_class,TSB_summary$pTNMStage_v8)[,2]),
                  table(TSB_summary$Zhao_class,TSB_summary$pTNMStage_v8)[,3]/sum(table(TSB_summary$Zhao_class,TSB_summary$pTNMStage_v8)[,3]))

#chi-squared test
chi.table <- TSB_summary
chi.table$Zhao_class <- as.character(chi.table$Zhao_class)
chi.table$Zhao_class[which(chi.table$Zhao_class %in% c("High","Low"))] <- "Non-discordant"
chisq_p <- chisq.test(table(chi.table$Zhao_class,chi.table$pTNMStage_v8))$p.value

ggplot(df,aes(Stage,risk.perc*100,fill=Risk)) + geom_col(width=0.8,position = "dodge") +
  ylab("Survival risk classification (%)")+
  xlab("")+ggtitle("Zhao")+
  geom_text(aes(y = risk.perc*100+1,label = paste0(signif(risk.perc*100,digits = 2),"%")),position = position_dodge(width=0.8),size=3)+
  geom_text(data=df, aes(x=2.5,y=75,label = paste0("chi-sq test, p=",signif(chisq_p,digits = 1))),size=5)+
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0, 0), limits = c(0,100))+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme(legend.position = "none",axis.text.x = element_text(size=12),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#Song risk-class frequency
df <- expand.grid( Risk = c("Low","High","Discordant"),Stage = c("Stage I","Stage II","Stage III"), perc = NA)
df$perc <- rep(table(TSB_summary$pTNMStage_v8)*100/nrow(TSB_summary),each=3)
df$risk.perc <- c(table(TSB_summary$Song_class,TSB_summary$pTNMStage_v8)[,1]/sum(table(TSB_summary$Song_class,TSB_summary$pTNMStage_v8)[,1]),
                  table(TSB_summary$Song_class,TSB_summary$pTNMStage_v8)[,2]/sum(table(TSB_summary$Song_class,TSB_summary$pTNMStage_v8)[,2]),
                  table(TSB_summary$Song_class,TSB_summary$pTNMStage_v8)[,3]/sum(table(TSB_summary$Song_class,TSB_summary$pTNMStage_v8)[,3]))

#chi-squared test
chi.table <- TSB_summary
chi.table$Song_class <- as.character(chi.table$Song_class)
chi.table$Song_class[which(chi.table$Song_class %in% c("High","Low"))] <- "Non-discordant"
chisq_p <- chisq.test(table(chi.table$Song_class,chi.table$pTNMStage_v8))$p.value

ggplot(df,aes(Stage,risk.perc*100,fill=Risk)) + geom_col(width=0.8,position = "dodge") +
  ylab("Survival risk classification (%)")+
  xlab("")+ggtitle("Song")+
  geom_text(aes(y = risk.perc*100+1,label = paste0(signif(risk.perc*100,digits = 2),"%")),position = position_dodge(width=0.8),size=3)+
  geom_text(data=df, aes(x=2.5,y=75,label = paste0("chi-sq test, p=",signif(chisq_p,digits = 1))),size=5)+
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0, 0), limits = c(0,100))+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme(legend.position = "none",axis.text.x = element_text(size=12),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#Jin risk-class frequency
df <- expand.grid( Risk = c("Low","High","Discordant"),Stage = c("Stage I","Stage II","Stage III"), perc = NA)
df$perc <- rep(table(TSB_summary$pTNMStage_v8)*100/nrow(TSB_summary),each=3)
df$risk.perc <- c(table(TSB_summary$Jin_class,TSB_summary$pTNMStage_v8)[,1]/sum(table(TSB_summary$Jin_class,TSB_summary$pTNMStage_v8)[,1]),
                  table(TSB_summary$Jin_class,TSB_summary$pTNMStage_v8)[,2]/sum(table(TSB_summary$Jin_class,TSB_summary$pTNMStage_v8)[,2]),
                  table(TSB_summary$Jin_class,TSB_summary$pTNMStage_v8)[,3]/sum(table(TSB_summary$Jin_class,TSB_summary$pTNMStage_v8)[,3]))

#chi-squared test
chi.table <- TSB_summary
chi.table$Jin_class <- as.character(chi.table$Jin_class)
chi.table$Jin_class[which(chi.table$Jin_class %in% c("High","Low"))] <- "Non-discordant"
chisq_p <- chisq.test(table(chi.table$Jin_class,chi.table$pTNMStage_v8))$p.value

ggplot(df,aes(Stage,risk.perc*100,fill=Risk)) + geom_col(width=0.8,position = "dodge") +
  ylab("Survival risk classification (%)")+
  xlab("")+ggtitle("Jin")+
  geom_text(aes(y = risk.perc*100+1,label = paste0(signif(risk.perc*100,digits = 2),"%")),position = position_dodge(width=0.8),size=3)+
  geom_text(data=df, aes(x=2.5,y=75,label = paste0("chi-sq test, p=",signif(chisq_p,digits = 1))),size=5)+
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0, 0), limits = c(0,100))+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme(legend.position = "none",axis.text.x = element_text(size=12),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#LiFeng risk-class frequency
df <- expand.grid( Risk = c("Low","High","Discordant"),Stage = c("Stage I","Stage II","Stage III"), perc = NA)
df$perc <- rep(table(TSB_summary$pTNMStage_v8)*100/nrow(TSB_summary),each=3)
df$risk.perc <- c(table(TSB_summary$LiFeng_class,TSB_summary$pTNMStage_v8)[,1]/sum(table(TSB_summary$LiFeng_class,TSB_summary$pTNMStage_v8)[,1]),
                  table(TSB_summary$LiFeng_class,TSB_summary$pTNMStage_v8)[,2]/sum(table(TSB_summary$LiFeng_class,TSB_summary$pTNMStage_v8)[,2]),
                  table(TSB_summary$LiFeng_class,TSB_summary$pTNMStage_v8)[,3]/sum(table(TSB_summary$LiFeng_class,TSB_summary$pTNMStage_v8)[,3]))

#chi-squared test
chi.table <- TSB_summary
chi.table$LiFeng_class <- as.character(chi.table$LiFeng_class)
chi.table$LiFeng_class[which(chi.table$LiFeng_class %in% c("High","Low"))] <- "Non-discordant"
chisq_p <- chisq.test(table(chi.table$LiFeng_class,chi.table$pTNMStage_v8))$p.value

ggplot(df,aes(Stage,risk.perc*100,fill=Risk)) + geom_col(width=0.8,position = "dodge") +
  ylab("Survival risk classification (%)")+
  xlab("")+ggtitle("LiFeng")+
  geom_text(aes(y = risk.perc*100+1,label = paste0(signif(risk.perc*100,digits = 2),"%")),position = position_dodge(width=0.8),size=3)+
  geom_text(data=df, aes(x=2.5,y=75,label = paste0("chi-sq test, p=",signif(chisq_p,digits = 1))),size=5)+
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0, 0), limits = c(0,100))+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme(legend.position = "none",axis.text.x = element_text(size=12),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))




### Extended Data Fig 3
clust_heatmap(mat2)
clust_heatmap(mat3)
clust_heatmap(mat4)
clust_heatmap(mat5)
clust_heatmap(mat6)
clust_heatmap(mat7)




### Extended Data Fig 4

## ED 4A
clust_heatmap(mat1)


## ED 4B

#number of patients
n_patient <- gsub("_.*",colnames(LUAD_vsd),replacement = "") %>% unique() %>% length()

#prepare df for concordant rate of each signature
sig_cluster <- data.frame(rate = rep(NA,n_patient*7),Signature=c(rep("ORACLE",n_patient),
                                                                 rep("Li JAMA Oncol 2017", n_patient),
                                                                 rep("Wang Front Immunol 2022", n_patient),
                                                                 rep("Zhao Lung Cancer 2020", n_patient),
                                                                 rep("Song Sci Rep 2022", n_patient),
                                                                 rep("Jin J Immunol Res 2022", n_patient),
                                                                 rep("Li Sci Rep 2022", n_patient)), cluster_no = rep(1:n_patient, 7))
#calculate concordant rate
for(i in 1:n_patient){
  sig_cluster$rate[i] <- clusterSigs(mat1,i)
  sig_cluster$rate[i+n_patient] <- clusterSigs(mat2,i)
  sig_cluster$rate[i+n_patient*2] <- clusterSigs(mat3,i)
  sig_cluster$rate[i+n_patient*3] <- clusterSigs(mat4,i)
  sig_cluster$rate[i+n_patient*4] <- clusterSigs(mat5,i)
  sig_cluster$rate[i+n_patient*5] <- clusterSigs(mat6,i)
  sig_cluster$rate[i+n_patient*6] <- clusterSigs(mat7,i)
}

#level signatures
sig_cluster$Signature <- factor(sig_cluster$Signature,levels = c("ORACLE","Li JAMA Oncol 2017","Wang Front Immunol 2022","Zhao Lung Cancer 2020","Song Sci Rep 2022","Jin J Immunol Res 2022","Li Sci Rep 2022"))

#line plot
ggplot(sig_cluster,aes(cluster_no,rate,col=Signature))+
  geom_line(size=1)+
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0,0),limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,120,20),expand = c(0,0),limits = c(0,n_patient))+
  scale_color_manual(values = c("#A6CEE3","#1F78B4","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00"))+
  geom_vline(xintercept = c(2,10,60),lty="dashed",col="gray")+
  xlab("No. of clusters") + ylab("Proportion of patients with\nall samples in the same cluster") + ggtitle(label = "Hclust assessment of tumour sampling bias")+
  theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))


## ED 4C

#expression sd
df1_sd <- ExpSD.df(mat1)
df2_sd <- ExpSD.df(mat2)
df3_sd <- ExpSD.df(mat3)
df4_sd <- ExpSD.df(mat4)
df5_sd <- ExpSD.df(mat5)
df6_sd <- ExpSD.df(mat6)
df7_sd <- ExpSD.df(mat7)

#mean sd
annot.df1 <- MeanSD.df(mat1,"ORACLE")
annot.df2 <- MeanSD.df(mat2,"Li JAMA Oncol 2017")
annot.df3 <- MeanSD.df(mat3,"Wang Front Immunol 2022")
annot.df4 <- MeanSD.df(mat4,"Zhao Lung Cancer 2020")
annot.df5 <- MeanSD.df(mat5,"Song Sci Rep 2022")
annot.df6 <- MeanSD.df(mat6,"Jin J Immunol Res 2022")
annot.df7 <- MeanSD.df(mat7,"Li Sci Rep 2022")

#combine results for all signatures
meta_mat <- cbind(df1_sd,df2_sd,df3_sd,df4_sd,df5_sd,df6_sd,df7_sd)
meta_annot <- rbind(annot.df1,annot.df2,annot.df3,annot.df4,annot.df5,annot.df6,annot.df7)

#colour
pal <- colorRampPalette(c(desaturate("steelblue4"), "steelblue4"))

#scale expression sd matrix
tmp <- as.data.frame(scale(t(meta_mat),center=T,scale=T))
rownames(meta_annot) <- rownames(tmp)

pheatmap(tmp,show_rownames = T , show_colnames = F,cluster_rows=F,cluster_cols=F,
         fontsize_col=10,border_color=NA,fontsize=4,legend=T,cellwidth=5,cellheight=2,
         color=wes_palette("Moonrise1", length(seq(floor(min(tmp)),ceiling(max(tmp)),by=0.5))-1, type="continuous"),breaks=seq(floor(min(tmp)),ceiling(max(tmp)),by=0.5),
         main='sd Expression',annotation_row = meta_annot,
         annotation_colors = list(signature = setNames(c("#A6CEE3","#1F78B4","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")
                                                       ,c("ORACLE",
                                                          "Li JAMA Oncol 2017",
                                                          "Wang Front Immunol 2022",
                                                          "Zhao Lung Cancer 2020",
                                                          "Song Sci Rep 2022",
                                                          "Jin J Immunol Res 2022",
                                                          "Li Sci Rep 2022")),
                                  mean_var = brewer.pal(9,'Purples')) 
)


## ED 4D

#retrieve mean of expression variation and Wilcoxon test compared with ORALCE
meta_annot$pval <- NA
meta_annot$pval[which(meta_annot$signature == "Jin J Immunol Res 2022")] <- wilcox.test(meta_annot$mean_var[which(meta_annot$signature == "ORACLE")],meta_annot$mean_var[which(meta_annot$signature == "Jin J Immunol Res 2022")],paired = F)$p.value
meta_annot$pval[which(meta_annot$signature == "Li JAMA Oncol 2017")] <- wilcox.test(meta_annot$mean_var[which(meta_annot$signature == "ORACLE")],meta_annot$mean_var[which(meta_annot$signature == "Li JAMA Oncol 2017")],paired = F)$p.value
meta_annot$pval[which(meta_annot$signature == "Wang Front Immunol 2022")] <- wilcox.test(meta_annot$mean_var[which(meta_annot$signature == "ORACLE")],meta_annot$mean_var[which(meta_annot$signature == "Wang Front Immunol 2022")],paired = F)$p.value
meta_annot$pval[which(meta_annot$signature == "Song Sci Rep 2022")] <- wilcox.test(meta_annot$mean_var[which(meta_annot$signature == "ORACLE")],meta_annot$mean_var[which(meta_annot$signature == "Song Sci Rep 2022")],paired = F)$p.value
meta_annot$pval[which(meta_annot$signature == "Li Sci Rep 2022")] <- wilcox.test(meta_annot$mean_var[which(meta_annot$signature == "ORACLE")],meta_annot$mean_var[which(meta_annot$signature == "Li Sci Rep 2022")],paired = F)$p.value
meta_annot$pval[which(meta_annot$signature == "Zhao Lung Cancer 2020")] <- wilcox.test(meta_annot$mean_var[which(meta_annot$signature == "ORACLE")],meta_annot$mean_var[which(meta_annot$signature == "Zhao Lung Cancer 2020")],paired = F)$p.value
meta_annot$signature <- factor(meta_annot$signature, levels = c("ORACLE","Jin J Immunol Res 2022","Li JAMA Oncol 2017","Wang Front Immunol 2022","Song Sci Rep 2022","Li Sci Rep 2022","Zhao Lung Cancer 2020"))

#assign p-value compare between signatures
pval.df <- meta_annot[which(meta_annot$signature != "ORACLE"),c(2:3)]
pval.df <- pval.df[which(!duplicated(pval.df$signature)),]
pval.df$group1 <- "ORACLE"
colnames(pval.df)[1] <- "group2"
pval.df$pval <- signif(pval.df$pval, digits = 2)

#order the signatures
pval.df$order <- NA
pval.df$order[which(pval.df$group2 == "Jin J Immunol Res 2022")] <- 1
pval.df$order[which(pval.df$group2 == "Li JAMA Oncol 2017")] <- 2
pval.df$order[which(pval.df$group2 == "Wang Front Immunol 2022")] <- 3
pval.df$order[which(pval.df$group2 == "Song Sci Rep 2022")] <- 4
pval.df$order[which(pval.df$group2 == "Shukla JNCI 2017")] <- 5
pval.df$order[which(pval.df$group2 == "Li Sci Rep 2022")] <- 6
pval.df$order[which(pval.df$group2 == "Zhao Lung Cancer 2020")] <- 7
pval.df <- arrange(pval.df, order)

#boxplot
ggplot(meta_annot)+geom_boxplot(aes(fct_reorder(signature,mean_var),mean_var,fill=signature),width=0.5, alpha = 0.8, outlier.alpha = 0)+
  geom_jitter(aes(fct_reorder(signature,mean_var),mean_var,fill=signature),pch = 21, col = "black", width = 0.2) +
  xlab("")+ylab("Mean variability")+
  scale_y_continuous(breaks = seq(0,4,1),expand = c(0,0),limits = c(0,4))+
  scale_fill_manual(values = c("#A6CEE3","#FDBF6F","#1F78B4","#33A02C","#E31A1C","#FF7F00","#FB9A99"))+
  stat_pvalue_manual(data = pval.df, y.position = 2.5, step.increase = 0.1, label = "pval") +
  theme(axis.text.x = element_text(angle = 90 ,size=5, vjust=0.5,hjust = 1),axis.ticks.x = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "none")


## ED 4E

#Illustration of Bachtiary method
exploratory_RiskScore <- Tx321_riskscore

#create df with ORACLE mean RS and ORACLE SD
exploratory_RiskScore$mean <- NA
for(i in 1:nrow(exploratory_RiskScore)){
  exploratory_RiskScore$mean[i] <- exploratory_RiskScore$ORACLE_riskscore[which(exploratory_RiskScore$cruk_id == exploratory_RiskScore$cruk_id[i])] %>% mean()
}

exploratory_RiskScore$sd <- NA
for(i in 1:nrow(exploratory_RiskScore)){
  exploratory_RiskScore$sd[i] <- exploratory_RiskScore$ORACLE_riskscore[which(exploratory_RiskScore$cruk_id == exploratory_RiskScore$cruk_id[i])] %>% sd()
}

exploratory_RiskScore$deviation <- NA
for(i in 1:nrow(exploratory_RiskScore)){
  exploratory_RiskScore$deviation[i] <- exploratory_RiskScore$ORACLE_riskscore[i] - exploratory_RiskScore$mean[i]
}

#keep only multi-region samples
exploratory_RiskScore <- exploratory_RiskScore[-which(is.na(exploratory_RiskScore$sd)),]
exploratory_RiskScore <- arrange(exploratory_RiskScore,sd)

#plot by difference from mean
ggplot(exploratory_RiskScore,aes(fct_reorder(cruk_id, sd, .fun=min),deviation))+
  geom_point(alpha = 0.6,pch=16,size = 4)+
  geom_line(col = "black")+
  geom_hline(yintercept = 0,linetype="dashed")+
  xlab("")+
  ylab("Risk-score\nDifference")+
  scale_y_continuous(breaks=seq(-1,1,0.5),expand=c(0,0),limits=c(-1,1))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))


## ED 4F

#calculate risk-score variations
var.df <- Tx321_riskscore[,c("cruk_id","ORACLE_riskscore","Li_RS","Wang_RS","Zhao_RS","Song_RS","Jin_RS","LiFeng_RS")]

for(i in 1:nrow(var.df)){
  var.df$ORACLE_var[i] <- Tx321_riskscore$ORACLE_riskscore[which(Tx321_riskscore$cruk_id == Tx321_riskscore$cruk_id[i])] %>% var()
  var.df$Li_var[i] <- Tx321_riskscore$Li_RS[which(Tx321_riskscore$cruk_id == Tx321_riskscore$cruk_id[i])] %>% var()
  var.df$Wang_var[i] <- Tx321_riskscore$Wang_RS[which(Tx321_riskscore$cruk_id == Tx321_riskscore$cruk_id[i])] %>% var()
  var.df$Zhao_var[i] <- Tx321_riskscore$Zhao_RS[which(Tx321_riskscore$cruk_id == Tx321_riskscore$cruk_id[i])] %>% var()
  var.df$Song_var[i] <- Tx321_riskscore$Song_RS[which(Tx321_riskscore$cruk_id == Tx321_riskscore$cruk_id[i])] %>% var()
  var.df$Jin_var[i] <- Tx321_riskscore$Jin_RS[which(Tx321_riskscore$cruk_id == Tx321_riskscore$cruk_id[i])] %>% var()
  var.df$LiFeng_var[i] <- Tx321_riskscore$LiFeng_RS[which(Tx321_riskscore$cruk_id == Tx321_riskscore$cruk_id[i])] %>% var()
}

#remove single-region samples
var.df <- var.df[which(!is.na(var.df$ORACLE_var)),]

#mean risk score
for(i in 1:nrow(var.df)){
  var.df$ORACLE_mean[i] <- var.df$ORACLE_riskscore[which(var.df$cruk_id == var.df$cruk_id[i])] %>% mean()
  var.df$Li_mean[i] <- var.df$Li_RS[which(var.df$cruk_id == var.df$cruk_id[i])] %>% mean()
  var.df$Wang_mean[i] <- var.df$Wang_RS[which(var.df$cruk_id == var.df$cruk_id[i])] %>% mean()
  var.df$Zhao_mean[i] <- var.df$Zhao_RS[which(var.df$cruk_id == var.df$cruk_id[i])] %>% mean()
  var.df$Song_mean[i] <- var.df$Song_RS[which(var.df$cruk_id == var.df$cruk_id[i])] %>% mean()
  var.df$Jin_mean[i] <- var.df$Jin_RS[which(var.df$cruk_id == var.df$cruk_id[i])] %>% mean()
  var.df$LiFeng_mean[i] <- var.df$LiFeng_RS[which(var.df$cruk_id == var.df$cruk_id[i])] %>% mean()
}

summarise.ratio.df <- var.df[which(!duplicated(var.df$cruk_id)),-c(2:8) ]

#DOI: 10.1158/1078-0432.CCR-06-0357
#"When more than one sample per patient is analyzed, the variance of the mean value per patient decreases as the number of replicates per patient increases."

#biop1
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_1[i] <- (mean(summarise.ratio.df$ORACLE_var)/1) / ((mean(summarise.ratio.df$ORACLE_var)/1) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_1[i] <- (mean(summarise.ratio.df$Li_var)/1) / ((mean(summarise.ratio.df$Li_var)/1) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_1[i] <- (mean(summarise.ratio.df$Wang_var)/1) / ((mean(summarise.ratio.df$Wang_var)/1) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_1[i] <- (mean(summarise.ratio.df$Zhao_var)/1) / ((mean(summarise.ratio.df$Zhao_var)/1) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_1[i] <- (mean(summarise.ratio.df$Jin_var)/1) / ((mean(summarise.ratio.df$Jin_var)/1) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_1[i] <- (mean(summarise.ratio.df$Song_var)/1) / ((mean(summarise.ratio.df$Song_var)/1) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_1[i] <- (mean(summarise.ratio.df$LiFeng_var)/1) / ((mean(summarise.ratio.df$LiFeng_var)/1) + var(summarise.ratio.df$LiFeng_mean))
}

#biop2
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_2[i] <- (mean(summarise.ratio.df$ORACLE_var)/2) / ((mean(summarise.ratio.df$ORACLE_var)/2) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_2[i] <- (mean(summarise.ratio.df$Li_var)/2) / ((mean(summarise.ratio.df$Li_var)/2) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_2[i] <- (mean(summarise.ratio.df$Wang_var)/2) / ((mean(summarise.ratio.df$Wang_var)/2) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_2[i] <- (mean(summarise.ratio.df$Zhao_var)/2) / ((mean(summarise.ratio.df$Zhao_var)/2) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_2[i] <- (mean(summarise.ratio.df$Jin_var)/2) / ((mean(summarise.ratio.df$Jin_var)/2) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_2[i] <- (mean(summarise.ratio.df$Song_var)/2) / ((mean(summarise.ratio.df$Song_var)/2) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_2[i] <- (mean(summarise.ratio.df$LiFeng_var)/2) / ((mean(summarise.ratio.df$LiFeng_var)/2) + var(summarise.ratio.df$LiFeng_mean))
}

#biop3
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_3[i] <- (mean(summarise.ratio.df$ORACLE_var)/3) / ((mean(summarise.ratio.df$ORACLE_var)/3) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_3[i] <- (mean(summarise.ratio.df$Li_var)/3) / ((mean(summarise.ratio.df$Li_var)/3) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_3[i] <- (mean(summarise.ratio.df$Wang_var)/3) / ((mean(summarise.ratio.df$Wang_var)/3) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_3[i] <- (mean(summarise.ratio.df$Zhao_var)/3) / ((mean(summarise.ratio.df$Zhao_var)/3) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_3[i] <- (mean(summarise.ratio.df$Jin_var)/3) / ((mean(summarise.ratio.df$Jin_var)/3) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_3[i] <- (mean(summarise.ratio.df$Song_var)/3) / ((mean(summarise.ratio.df$Song_var)/3) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_3[i] <- (mean(summarise.ratio.df$LiFeng_var)/3) / ((mean(summarise.ratio.df$LiFeng_var)/3) + var(summarise.ratio.df$LiFeng_mean))
}

#biop4
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_4[i] <- (mean(summarise.ratio.df$ORACLE_var)/4) / ((mean(summarise.ratio.df$ORACLE_var)/4) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_4[i] <- (mean(summarise.ratio.df$Li_var)/4) / ((mean(summarise.ratio.df$Li_var)/4) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_4[i] <- (mean(summarise.ratio.df$Wang_var)/4) / ((mean(summarise.ratio.df$Wang_var)/4) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_4[i] <- (mean(summarise.ratio.df$Zhao_var)/4) / ((mean(summarise.ratio.df$Zhao_var)/4) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_4[i] <- (mean(summarise.ratio.df$Jin_var)/4) / ((mean(summarise.ratio.df$Jin_var)/4) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_4[i] <- (mean(summarise.ratio.df$Song_var)/4) / ((mean(summarise.ratio.df$Song_var)/4) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_4[i] <- (mean(summarise.ratio.df$LiFeng_var)/4) / ((mean(summarise.ratio.df$LiFeng_var)/4) + var(summarise.ratio.df$LiFeng_mean))
}

#biop5
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_5[i] <- (mean(summarise.ratio.df$ORACLE_var)/5) / ((mean(summarise.ratio.df$ORACLE_var)/5) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_5[i] <- (mean(summarise.ratio.df$Li_var)/5) / ((mean(summarise.ratio.df$Li_var)/5) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_5[i] <- (mean(summarise.ratio.df$Wang_var)/5) / ((mean(summarise.ratio.df$Wang_var)/5) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_5[i] <- (mean(summarise.ratio.df$Zhao_var)/5) / ((mean(summarise.ratio.df$Zhao_var)/5) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_5[i] <- (mean(summarise.ratio.df$Jin_var)/5) / ((mean(summarise.ratio.df$Jin_var)/5) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_5[i] <- (mean(summarise.ratio.df$Song_var)/5) / ((mean(summarise.ratio.df$Song_var)/5) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_5[i] <- (mean(summarise.ratio.df$LiFeng_var)/5) / ((mean(summarise.ratio.df$LiFeng_var)/5) + var(summarise.ratio.df$LiFeng_mean))
}

#biop6
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_6[i] <- (mean(summarise.ratio.df$ORACLE_var)/6) / ((mean(summarise.ratio.df$ORACLE_var)/6) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_6[i] <- (mean(summarise.ratio.df$Li_var)/6) / ((mean(summarise.ratio.df$Li_var)/6) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_6[i] <- (mean(summarise.ratio.df$Wang_var)/6) / ((mean(summarise.ratio.df$Wang_var)/6) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_6[i] <- (mean(summarise.ratio.df$Zhao_var)/6) / ((mean(summarise.ratio.df$Zhao_var)/6) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_6[i] <- (mean(summarise.ratio.df$Jin_var)/6) / ((mean(summarise.ratio.df$Jin_var)/6) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_6[i] <- (mean(summarise.ratio.df$Song_var)/6) / ((mean(summarise.ratio.df$Song_var)/6) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_6[i] <- (mean(summarise.ratio.df$LiFeng_var)/6) / ((mean(summarise.ratio.df$LiFeng_var)/6) + var(summarise.ratio.df$LiFeng_mean))
}

#biop7
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_7[i] <- (mean(summarise.ratio.df$ORACLE_var)/7) / ((mean(summarise.ratio.df$ORACLE_var)/7) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_7[i] <- (mean(summarise.ratio.df$Li_var)/7) / ((mean(summarise.ratio.df$Li_var)/7) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_7[i] <- (mean(summarise.ratio.df$Wang_var)/7) / ((mean(summarise.ratio.df$Wang_var)/7) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_7[i] <- (mean(summarise.ratio.df$Zhao_var)/7) / ((mean(summarise.ratio.df$Zhao_var)/7) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_7[i] <- (mean(summarise.ratio.df$Jin_var)/7) / ((mean(summarise.ratio.df$Jin_var)/7) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_7[i] <- (mean(summarise.ratio.df$Song_var)/7) / ((mean(summarise.ratio.df$Song_var)/7) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_7[i] <- (mean(summarise.ratio.df$LiFeng_var)/7) / ((mean(summarise.ratio.df$LiFeng_var)/7) + var(summarise.ratio.df$LiFeng_mean))
}

#biop8
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_8[i] <- (mean(summarise.ratio.df$ORACLE_var)/8) / ((mean(summarise.ratio.df$ORACLE_var)/8) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_8[i] <- (mean(summarise.ratio.df$Li_var)/8) / ((mean(summarise.ratio.df$Li_var)/8) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_8[i] <- (mean(summarise.ratio.df$Wang_var)/8) / ((mean(summarise.ratio.df$Wang_var)/8) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_8[i] <- (mean(summarise.ratio.df$Zhao_var)/8) / ((mean(summarise.ratio.df$Zhao_var)/8) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_8[i] <- (mean(summarise.ratio.df$Jin_var)/8) / ((mean(summarise.ratio.df$Jin_var)/8) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_8[i] <- (mean(summarise.ratio.df$Song_var)/8) / ((mean(summarise.ratio.df$Song_var)/8) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_8[i] <- (mean(summarise.ratio.df$LiFeng_var)/8) / ((mean(summarise.ratio.df$LiFeng_var)/8) + var(summarise.ratio.df$LiFeng_mean))
}

#biop9
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_9[i] <- (mean(summarise.ratio.df$ORACLE_var)/9) / ((mean(summarise.ratio.df$ORACLE_var)/9) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_9[i] <- (mean(summarise.ratio.df$Li_var)/9) / ((mean(summarise.ratio.df$Li_var)/9) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_9[i] <- (mean(summarise.ratio.df$Wang_var)/9) / ((mean(summarise.ratio.df$Wang_var)/9) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_9[i] <- (mean(summarise.ratio.df$Zhao_var)/9) / ((mean(summarise.ratio.df$Zhao_var)/9) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_9[i] <- (mean(summarise.ratio.df$Jin_var)/9) / ((mean(summarise.ratio.df$Jin_var)/9) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_9[i] <- (mean(summarise.ratio.df$Song_var)/9) / ((mean(summarise.ratio.df$Song_var)/9) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_9[i] <- (mean(summarise.ratio.df$LiFeng_var)/9) / ((mean(summarise.ratio.df$LiFeng_var)/9) + var(summarise.ratio.df$LiFeng_mean))
}

#biop10
for(i in 1:nrow(summarise.ratio.df)){
  summarise.ratio.df$ORACLE_WT_10[i] <- (mean(summarise.ratio.df$ORACLE_var)/10) / ((mean(summarise.ratio.df$ORACLE_var)/10) + var(summarise.ratio.df$ORACLE_mean))
  summarise.ratio.df$Li_WT_10[i] <- (mean(summarise.ratio.df$Li_var)/10) / ((mean(summarise.ratio.df$Li_var)/10) + var(summarise.ratio.df$Li_mean))
  summarise.ratio.df$Wang_WT_10[i] <- (mean(summarise.ratio.df$Wang_var)/10) / ((mean(summarise.ratio.df$Wang_var)/10) + var(summarise.ratio.df$Wang_mean))
  summarise.ratio.df$Zhao_WT_10[i] <- (mean(summarise.ratio.df$Zhao_var)/10) / ((mean(summarise.ratio.df$Zhao_var)/10) + var(summarise.ratio.df$Zhao_mean))
  summarise.ratio.df$Jin_WT_10[i] <- (mean(summarise.ratio.df$Jin_var)/10) / ((mean(summarise.ratio.df$Jin_var)/10) + var(summarise.ratio.df$Jin_mean))
  summarise.ratio.df$Song_WT_10[i] <- (mean(summarise.ratio.df$Song_var)/10) / ((mean(summarise.ratio.df$Song_var)/10) + var(summarise.ratio.df$Song_mean))
  summarise.ratio.df$LiFeng_WT_10[i] <- (mean(summarise.ratio.df$LiFeng_var)/10) / ((mean(summarise.ratio.df$LiFeng_var)/10) + var(summarise.ratio.df$LiFeng_mean))
}

#line plot
ggplot(biopsy_need_df, aes(factor(biopsy), value)) + geom_point(aes(col = signature,pch=label),size=3) +
  ylab("W/T per #biopsy") + xlab("#biopsy") +
  geom_line(aes(factor(biopsy), value, group=signature,col=signature)) +
  scale_color_manual(values = c("#A6CEE3","#1F78B4","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")) + 
  geom_hline(yintercept = 0.15, lty="dashed",col="darkblue") +
  scale_shape_manual(values = c(21,19)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.4),position = "right") +
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))

#least biopsy needed for each signature using linear fit between points of biopsy 1 and 2
md <- lm(biopsy_need_df$value[which(biopsy_need_df$signature=="ORACLE" & biopsy_need_df$biopsy%in%c(1,2))]~c(1,2))
r <- summary(md)
biop_ORACLE <- (0.15-r$coefficients[1,1])/r$coefficients[2,1]

md <- lm(biopsy_need_df$value[which(biopsy_need_df$signature=="Li JAMA Oncol 2017" & biopsy_need_df$biopsy%in%c(2,3))]~c(2,3))
r <- summary(md)
biop_Li <- (0.15-r$coefficients[1,1])/r$coefficients[2,1]

md <- lm(biopsy_need_df$value[which(biopsy_need_df$signature=="Wang Front Immunol 2022" & biopsy_need_df$biopsy%in%c(1,2))]~c(1,2))
r <- summary(md)
biop_Wang <- (0.15-r$coefficients[1,1])/r$coefficients[2,1]

md <- lm(biopsy_need_df$value[which(biopsy_need_df$signature=="Zhao Lung Cancer 2020" & biopsy_need_df$biopsy%in%c(2,3))]~c(2,3))
r <- summary(md)
biop_Zhao <- (0.15-r$coefficients[1,1])/r$coefficients[2,1]

md <- lm(biopsy_need_df$value[which(biopsy_need_df$signature=="Song Sci Rep 2022" & biopsy_need_df$biopsy%in%c(1,2))]~c(1,2))
r <- summary(md)
biop_Song <- (0.15-r$coefficients[1,1])/r$coefficients[2,1]

md <- lm(biopsy_need_df$value[which(biopsy_need_df$signature=="Jin J Immunol Res 2022" & biopsy_need_df$biopsy%in%c(2,3))]~c(2,3))
r <- summary(md)
biop_Jin <- (0.15-r$coefficients[1,1])/r$coefficients[2,1]

md <- lm(biopsy_need_df$value[which(biopsy_need_df$signature=="Li Sci Rep 2022" & biopsy_need_df$biopsy%in%c(1,2))]~c(1,2))
r <- summary(md)
biop_LiFeng <- (0.15-r$coefficients[1,1])/r$coefficients[2,1]

least.biopsy <- data.frame(signature = c("ORACLE","Li JAMA Oncol 2017","Wang Front Immunol 2022","Zhao Lung Cancer 2020","Song Sci Rep 2022","Jin J Immunol Res 2022","Li Sci Rep 2022"),
                           biop_need = c(biop_ORACLE,biop_Li,biop_Wang,biop_Zhao,biop_Song,biop_Jin,biop_LiFeng))




### Extended Data Fig 5

## ED 5A

MVA.cox.res <- coxph(formula = Surv(Survdata_stageI$os_time, Survdata_stageI$cens_os) ~ 
                       Survdata_stageI$ORACLE_mean + 
                       Survdata_stageI$sex +
                       Survdata_stageI$age + 
                       Survdata_stageI$pack_years_calculated +  
                       Survdata_stageI$adjuvant_treatment_YN +
                       Survdata_stageI$TNM_combine +
                       Survdata_stageI$IASLC_LUAD_grade
)

MVA.cox.res <- summary(MVA.cox.res)
MVA.cox.res <- data.frame(HR = MVA.cox.res$coefficients[,2], lower_ci = MVA.cox.res$conf.int[,3], upper_ci = MVA.cox.res$conf.int[,4], P = MVA.cox.res$coefficients[,5])
MVA.cox.res$Predictors <- c("Mean RS","Female","Age","pack_years","Adjuvant treatment","Stage II","Stage III","Grade 2","Grade 3")

#add ref
idx <- min(grep("Female", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","Male")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Adjuvant", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","No adjuvant")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Stage", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","Stage I")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Grade", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","IMA")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

#turn data from chr to numeric
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#tidy up presenting values - round up and 95% CI
MVA.cox.res$pvalue <-signif(MVA.cox.res$P,digits = 1) 

#add haeder
MVA.cox.res <- rbind(c(1,NA,NA,"","","p-value","HR","CI"),MVA.cox.res)
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#Forest plot
forestplot(MVA.cox.res[,c(5:6)],graph.pos=2,
           mean = MVA.cox.res$HR,
           lower = MVA.cox.res$lower_ci,
           upper = MVA.cox.res$upper_ci,
           xlog=TRUE,title="Hazard Ratio",
           boxsize = 0.4,ci.vertices=TRUE,ci.vertices.height = 0.2,
           txt_gp=fpTxtGp(label=gpar(cex=0.9),
                          ticks=gpar(cex=0.9),
                          xlab=gpar(cex = 1),
                          title=gpar(cex = 1.2)),lwd.ci=1,colgap=unit(6,"mm"),
           col=fpColors(box="black", lines="black", zero = "gray50"),
           graphwidth = unit(12,"cm"))


## ED 5B

Sankey.df <- as.data.frame(table(Survdata_stageI$pTNMStage_v8,Survdata_stageI$ORACLE_class))
colnames(Sankey.df) <- c("Stage","Risk","Freq")
Sankey.df$Risk <- factor(Sankey.df$Risk,levels = c("Low","Discordant","High"))

ggplot(as.data.frame(Sankey.df), aes(y = Freq, axis1 = Stage, axis2 = Risk)) +
  geom_alluvium(aes(fill = Risk), width = 1/4) +
  scale_linetype_manual(values = c("blank", "solid"))+
  geom_stratum(width = 1/4,fill = "white")+
  geom_text(stat = "stratum",aes(label = stat(stratum)),vjust = 0.5,size = 5,direction = "y")+
  scale_x_discrete(limits = c("Stage", "ORACLE"), expand = c(0.05,0.05))+
  scale_fill_manual(values = c("#3B4992FF","azure4","#EE0000FF"))+
  scale_y_continuous(labels = NULL,expand = c(0,0))+
  ylab(NULL)+
  theme(legend.position = "none",legend.title = element_text(size=14),legend.text = element_text(size = 12), axis.text.x = element_text(size = 12),axis.ticks = element_blank(),panel.background = element_blank())


## ED 5C

#mircoarray expression matrices
ma_shedden <- column_to_rownames(ma_shedden,var = "PatientID")
ma_rousseaux <- column_to_rownames(ma_rousseaux,var = "PatientID")
ma_okayama <- column_to_rownames(ma_okayama,var = "PatientID")
ma_der <- column_to_rownames(ma_der,var = "PatientID")

ma_shedden <- as.data.frame(t(ma_shedden))
ma_rousseaux <- as.data.frame(t(ma_rousseaux))
ma_okayama <- as.data.frame(t(ma_okayama))
ma_der <- as.data.frame(t(ma_der))

#select oracle gene
ma_shedden_oracle <- ma_shedden[which(rownames(ma_shedden) %in% Supp_Table_5_ORACLE$Gene.Symbol), ]
ma_rousseaux_oracle <- ma_rousseaux[which(rownames(ma_rousseaux) %in% Supp_Table_5_ORACLE$Gene.Symbol), ]
ma_okayama_oracle <- ma_okayama[which(rownames(ma_okayama) %in% Supp_Table_5_ORACLE$Gene.Symbol), ]
ma_der_oracle <- ma_der[which(rownames(ma_der) %in% Supp_Table_5_ORACLE$Gene.Symbol), ]

#calculate ORACLE score
ORACLE <- Supp_Table_5_ORACLE[which(Supp_Table_5_ORACLE$Gene.Symbol %in% rownames(ma_shedden_oracle)),]
shedden_oracle_rs <- data.frame(RiskScore = colSums(ma_shedden_oracle[ORACLE$Gene.Symbol,] * ORACLE$Model.Coefficient), PatientID = colnames(ma_shedden_oracle))
shedden_oracle_rs <- left_join(shedden_oracle_rs, surv_shedden, by = "PatientID")

ORACLE <- Supp_Table_5_ORACLE[which(Supp_Table_5_ORACLE$Gene.Symbol %in% rownames(ma_rousseaux_oracle)),]
rousseaux_oracle_rs <- data.frame(RiskScore = colSums(ma_rousseaux_oracle[ORACLE$Gene.Symbol,] * ORACLE$Model.Coefficient), PatientID = colnames(ma_rousseaux_oracle))
rousseaux_oracle_rs <- left_join(rousseaux_oracle_rs, surv_rousseaux, by = "PatientID")

ORACLE <- Supp_Table_5_ORACLE[which(Supp_Table_5_ORACLE$Gene.Symbol %in% rownames(ma_okayama_oracle)),]
okayama_oracle_rs <- data.frame(RiskScore = colSums(ma_okayama_oracle[ORACLE$Gene.Symbol,] * ORACLE$Model.Coefficient), PatientID = colnames(ma_okayama_oracle))
okayama_oracle_rs <- left_join(okayama_oracle_rs, surv_okayama, by = "PatientID")

ORACLE <- Supp_Table_5_ORACLE[which(Supp_Table_5_ORACLE$Gene.Symbol %in% rownames(ma_der_oracle)),]
der_oracle_rs <- data.frame(RiskScore = colSums(ma_der_oracle[ORACLE$Gene.Symbol,] * ORACLE$Model.Coefficient), PatientID = colnames(ma_der_oracle))
der_oracle_rs <- left_join(der_oracle_rs, surv_der, by = "PatientID")

#UVA in stage I
coxmd <- coxph(Surv(time, status) ~ RiskScore , data = shedden_oracle_rs[which(shedden_oracle_rs$Histology == "LUAD" & shedden_oracle_rs$Stage %in% c("1a","1b")),])
Shedden_cox <- summary(coxmd)
Shedden_cox <- data.frame(HR=Shedden_cox$coefficients[,2], lower_ci=Shedden_cox$conf.int[,3], upper_ci=Shedden_cox$conf.int[,4], P=Shedden_cox$coefficients[,5], COEF=Shedden_cox$coefficients[1], seCOEF=Shedden_cox$coefficients[3])
Shedden_cox$cohort <- "Shedden"

coxmd <- coxph(Surv(time, status) ~ RiskScore , data = rousseaux_oracle_rs[which(rousseaux_oracle_rs$Histology == "LUAD" & rousseaux_oracle_rs$Stage %in% c("1a","1b")),])
Rousseaux_cox <- summary(coxmd)
Rousseaux_cox <- data.frame(HR=Rousseaux_cox$coefficients[,2], lower_ci=Rousseaux_cox$conf.int[,3], upper_ci=Rousseaux_cox$conf.int[,4], P=Rousseaux_cox$coefficients[,5], COEF=Rousseaux_cox$coefficients[1], seCOEF=Rousseaux_cox$coefficients[3])
Rousseaux_cox$cohort <- "Rousseaux"

coxmd <- coxph(Surv(time, status) ~ RiskScore , data = okayama_oracle_rs[which(okayama_oracle_rs$Histology == "LUAD" & okayama_oracle_rs$Stage %in% c("1A","1B")),])
Okayama_cox <- summary(coxmd)
Okayama_cox <- data.frame(HR=Okayama_cox$coefficients[,2], lower_ci=Okayama_cox$conf.int[,3], upper_ci=Okayama_cox$conf.int[,4], P=Okayama_cox$coefficients[,5], COEF=Okayama_cox$coefficients[1], seCOEF=Okayama_cox$coefficients[3])
Okayama_cox$cohort <- "Okayama"

coxmd <- coxph(Surv(time, status) ~ RiskScore , data = der_oracle_rs[which(der_oracle_rs$Histology == "LUAD" & der_oracle_rs$Stage %in% c("1a","1b")),])
Der_cox <- summary(coxmd)
Der_cox <- data.frame(HR=Der_cox$coefficients[,2], lower_ci=Der_cox$conf.int[,3], upper_ci=Der_cox$conf.int[,4], P=Der_cox$coefficients[,5], COEF=Der_cox$coefficients[1], seCOEF=Der_cox$coefficients[3])
Der_cox$cohort <- "Der"

#weighted HR
combine.cox <- rbind(Der_cox,Okayama_cox,Rousseaux_cox,Shedden_cox)

{
  # input: Cox regression hazard ratios and P-values for Uppsala + microarray cohorts
  tmp <- combine.cox
  rownames(tmp) <- tmp$cohort
  
  # use log values 
  tmp <- tmp[,5:6]
  sum.estimate <- meta.summaries(tmp[,1],tmp[,2],logscale=T)
  
  # meta-analysis plot
  metaplot(tmp[,1],tmp[,2], labels=rownames(tmp), xlab='Log Hazard Ratio', main='ORACLE', summn = sum.estimate$summary, sumse = sum.estimate$se.summary, sumnn= 1/sum.estimate$se.summary^2, xlim=c(-1, 5), zero=0, colors=meta.colors(box="#3182bd",lines="#a50f15", zero="red", summary="black",text="black"),xaxt='n')
  axis(1, at=log(c(0.5,1,2,4,8,16,32,64)), labels=c(0.5,1,2,4,8,16,32,64))
  
  # p-value, : https://www.bmj.com/content/343/bmj.d2304
  #If the upper and lower limits of a 95% CI are u and l respectively:
  #1 calculate the standard error: SE = (u − l)/(2×1.96)
  # 2 calculate the test statistic: z = Est/SE
  # 3 calculate the P value: P = exp(−0.717×z − 0.416×z2).
  Est <- 3.43 #sum.estimate$summary
  l <- 2.19
  u <- 5.37
  SE <- (u-l)/(2*1.96)
  z = Est/SE
  P = exp(-0.717*z - 0.416*z^2)
  # P= 2.841393e-05
}


## ED 5D

MVA.cox.res <- coxph(formula = Surv(Tx321_tumour_Survdata$lung_specific_time, Tx321_tumour_Survdata$cens_lung_specific) ~ 
                       Tx321_tumour_Survdata$ORACLE_mean + 
                       Tx321_tumour_Survdata$sex +
                       Tx321_tumour_Survdata$age + 
                       Tx321_tumour_Survdata$pack_years_calculated +  
                       Tx321_tumour_Survdata$adjuvant_treatment_YN +
                       Tx321_tumour_Survdata$TNM_combine 
)

MVA.cox.res <- summary(MVA.cox.res)
MVA.cox.res <- data.frame(HR = MVA.cox.res$coefficients[,2], lower_ci = MVA.cox.res$conf.int[,3], upper_ci = MVA.cox.res$conf.int[,4], P = MVA.cox.res$coefficients[,5])
MVA.cox.res$Predictors <- c("Mean RS","Female","Age","pack_years","Adjuvant treatment","Stage II","Stage III")

#add ref
idx <- min(grep("Female", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","Male")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Adjuvant", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","No adjuvant")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Stage", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","Stage I")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

#turn data from chr to numeric
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#tidy up presenting values - round up and 95% CI
MVA.cox.res$pvalue <-signif(MVA.cox.res$P,digits = 1) 

#add haeder
MVA.cox.res <- rbind(c(1,NA,NA,"","","p-value","HR","CI"),MVA.cox.res)
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#Forest plot
forestplot(MVA.cox.res[,c(5:6)],graph.pos=2,
           mean = MVA.cox.res$HR,
           lower = MVA.cox.res$lower_ci,
           upper = MVA.cox.res$upper_ci,
           xlog=TRUE,title="Hazard Ratio",
           boxsize = 0.4,ci.vertices=TRUE,ci.vertices.height = 0.2,
           txt_gp=fpTxtGp(label=gpar(cex=0.9),
                          ticks=gpar(cex=0.9),
                          xlab=gpar(cex = 1),
                          title=gpar(cex = 1.2)),lwd.ci=1,colgap=unit(6,"mm"),
           col=fpColors(box="black", lines="black", zero = "gray50"),
           graphwidth = unit(12,"cm"))


## ED 5E

MVA.cox.res <- coxph(formula = Surv(Relapse_survdata$os_time, Relapse_survdata$cens_os) ~ 
                       Relapse_survdata$ORACLE_mean + 
                       Relapse_survdata$sex +
                       Relapse_survdata$age + 
                       Relapse_survdata$pack_years_calculated +  
                       Relapse_survdata$adjuvant_treatment_YN +
                       Relapse_survdata$TNM_combine 
)

MVA.cox.res <- summary(MVA.cox.res)
MVA.cox.res <- data.frame(HR = MVA.cox.res$coefficients[,2], lower_ci = MVA.cox.res$conf.int[,3], upper_ci = MVA.cox.res$conf.int[,4], P = MVA.cox.res$coefficients[,5])
MVA.cox.res$Predictors <- c("Mean RS","Female","Age","pack_years","Adjuvant treatment","Stage II","Stage III")

#add ref
idx <- min(grep("Female", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","Male")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Adjuvant", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","No adjuvant")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

idx <- min(grep("Stage", MVA.cox.res$Predictors)) - 1
newRow <- c(1, 1, 1, "","Stage I")
MVA.cox.res <- rbind(MVA.cox.res[1:idx, ], newRow, MVA.cox.res[-(1:idx), ])

#turn data from chr to numeric
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#tidy up presenting values - round up and 95% CI
MVA.cox.res$pvalue <-signif(MVA.cox.res$P,digits = 1) 

#add haeder
MVA.cox.res <- rbind(c(1,NA,NA,"","","p-value","HR","CI"),MVA.cox.res)
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#Forest plot
forestplot(MVA.cox.res[,c(5:6)],graph.pos=2,
           mean = MVA.cox.res$HR,
           lower = MVA.cox.res$lower_ci,
           upper = MVA.cox.res$upper_ci,
           xlog=TRUE,title="Hazard Ratio",
           boxsize = 0.4,ci.vertices=TRUE,ci.vertices.height = 0.2,
           txt_gp=fpTxtGp(label=gpar(cex=0.9),
                          ticks=gpar(cex=0.9),
                          xlab=gpar(cex = 1),
                          title=gpar(cex = 1.2)),lwd.ci=1,colgap=unit(6,"mm"),
           col=fpColors(box="black", lines="black", zero = "gray50"),
           graphwidth = unit(12,"cm"))




### Extended Data Fig 6

## 6B

tmp <- IC50_results[which(IC50_results$significance == "TRUE"),]
sig.drugs <- fct_reorder(rownames(tmp), -tmp$estimate)
sig.drugs <- levels(sig.drugs)

for(i in 1:length(sig.drugs)){
  md <- cor.test(LUAD_IC50[,sig.drugs[i]],LUAD_IC50$riskscore,method="spearman",exact=FALSE)
  volc_p <- ggplot(LUAD_IC50, aes(LUAD_IC50[,sig.drugs[i]], riskscore)) + geom_point(size=5) + geom_smooth(method="lm",se=F) +
    xlab("Log(IC50)") + ylab("Risk score") + ggtitle(sig.drugs[i]) +
    scale_y_continuous(expand = c(0,0),limits = c(10.5,13.5)) +
    geom_text(aes(x=min(LUAD_IC50[,sig.drugs[i]], na.rm=T) + 1.4, y=13.4, label = paste0("R=",signif(md$estimate,digits = 2),", P=",signif(md$p.value,digits = 2)))) +
    theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())
  
  print(volc_p)
}




### Extended Data Fig 7

Tx321_tumour_Survdata <- left_join(Tx321_tumour_Survdata, all_tumour_df_20220209[,c("tumour_id_mphase_cruk","pN_stage")], by = "tumour_id_mphase_cruk")
Survdata_node_neg <- Tx321_tumour_Survdata[which(Tx321_tumour_Survdata$pN_stage == "0" | Tx321_tumour_Survdata$pN_stage == "X"),]
Survdata_node_pos <- Tx321_tumour_Survdata[which(Tx321_tumour_Survdata$pN_stage != "0" & Tx321_tumour_Survdata$pN_stage != "X"),]

#Node negative
Survdata_node_neg_noAD <- Survdata_node_neg[which(Survdata_node_neg$adjuvant_treatment_YN == "No adjuvant"),]
Survdata_node_neg_noAD$ORACLE_class <- factor(Survdata_node_neg_noAD$ORACLE_class, levels = c("Low","Discordant","High"))
j <- ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Survdata_node_neg_noAD),
                legend.labs = c("Low", "Discor","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                title="No adjuvant therapy", 
                legend="none", xlab="Time (Years)", ylab="Overall Survival",
                pval = T,pval.size = 8,
                ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())


j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j

Survdata_node_neg_withAD <- Survdata_node_neg[which(Survdata_node_neg$adjuvant_treatment_YN == "Adjuvant"),]
Survdata_node_neg_withAD$ORACLE_class <- factor(Survdata_node_neg_withAD$ORACLE_class, levels = c("Low","Discordant","High"))
j <- ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Survdata_node_neg_withAD),
                legend.labs = c("Low", "Discor","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                title="With adjuvant therapy", 
                legend="none", xlab="Time (Years)", ylab="Overall Survival",
                pval = T,pval.size = 8,
                ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j

#Node positive
Survdata_node_pos_noAD <- Survdata_node_pos[which(Survdata_node_pos$adjuvant_treatment_YN == "No adjuvant"),]
Survdata_node_pos_noAD$ORACLE_class <- factor(Survdata_node_pos_noAD$ORACLE_class, levels = c("Low","Discordant","High"))
j <- ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Survdata_node_pos_noAD),
                legend.labs = c("Low", "Discor","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                title="No adjuvant therapy", 
                legend="none", xlab="Time (Years)", ylab="Overall Survival",
                pval = T,pval.size = 8,
                ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j

Survdata_node_pos_withAD <- Survdata_node_pos[which(Survdata_node_pos$adjuvant_treatment_YN == "Adjuvant"),]
Survdata_node_pos_withAD$ORACLE_class <- factor(Survdata_node_pos_withAD$ORACLE_class, levels = c("Low","Discordant","High"))
j <- ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Survdata_node_pos_withAD),
                legend.labs = c("Low", "Discor","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                title="With adjuvant therapy", 
                legend="none", xlab="Time (Years)", ylab="Overall Survival",
                pval = T,pval.size = 8,
                ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())



j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j




### Extended Data Fig 8

ggplot(Features.df, aes(age, mean)) + geom_point() + xlab("Age") + ylab("ORACLE score") +
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(sex, mean)) + geom_boxplot() + xlab("Sex") + ylab("ORACLE score") + 
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(smoking_status, mean)) + geom_boxplot() + xlab("Smoking status") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(pTNMStage_v8, mean)) + geom_boxplot() + xlab("TNM stage") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(factor(biopsy_no), mean)) + geom_boxplot() + xlab("#Biopsy") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(Volume, mean)) + geom_point() + xlab("Volume") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(KI67_score, mean)) + geom_point() + xlab("Ki67") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(frac_abberant_genom_subcl, mean)) + geom_point() + xlab("SCNA-ITH") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(wFLOH, mean)) + geom_point() + xlab("FLOH") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(num_clonal_gds, mean)) + geom_point() + xlab("Clonal WGD") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(num_subclonal_gds, mean)) + geom_point() + xlab("Subclonal WGD") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(num_clonal_drivers, mean)) + geom_point() + xlab("Clonal drivers") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(num_subclonal_drivers, mean)) + geom_point() + xlab("Subclonal drivers") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())

ggplot(Features.df, aes(RecentSubclExpansionScore, mean)) + geom_point() + xlab("Recent Subclonal Expansion") + ylab("ORACLE score") +  
  scale_y_continuous(expand = c(0,0),limits = c(8,12)) +
  theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())




### Extended Data Fig 9

## ED 9A

#retreive clonal mutation
Clonal.SNV <- Tx421_regionMut[which(Tx421_regionMut$PyCloneClonal_SC == "C"),c("RegionID","PyCloneClonal_SC","Hugo_Symbol")]
Clonal.SNV <- as.data.frame(dcast(Clonal.SNV, RegionID~Hugo_Symbol, value.var = "PyCloneClonal_SC"))
colnames(Clonal.SNV)[1] <- "sample_name"

#binary df for all primary Tx samples, 0=no mutation, 1=mutation
clonal.freq.df <- data.frame(sample_name = gsub("LTX0",Tx421_riskscore$sample_name,replacement = "LTX"), bin = Tx421_riskscore$bin)
clonal.freq.df <- left_join(clonal.freq.df, Clonal.SNV)
for(i in 3:ncol(clonal.freq.df)){
  clonal.freq.df[which(is.na(clonal.freq.df[,i])),i] <- 0
}

for(i in 3:ncol(clonal.freq.df)){
  clonal.freq.df[which(clonal.freq.df[,i] > 1),i] <- 1
}

#remove low event number sample
rownames(clonal.freq.df) <- clonal.freq.df$sample_name
clonal.freq.df <- dplyr::select(clonal.freq.df, -contains("sample_name"))
clonal.freq.df <- clonal.freq.df[, -c(which(colSums(clonal.freq.df[,-1]) == 0) + 1 )]
clonal.freq.df <- clonal.freq.df[,-c(which(colSums(clonal.freq.df[,-1]) < nrow(clonal.freq.df)*0.02)+1)]
clonal.freq.df$bin <- factor(clonal.freq.df$bin, levels = c("Low","High"))

#Fisher exact test pvalue
clonal.result.df <- data.frame(Driver = colnames(clonal.freq.df)[-1],pvalue=NA,odds_ratio = NA)

for(i in 1:nrow(clonal.result.df)){
  clonal.result.df$pvalue[i] <- fisher.test(table(clonal.freq.df$bin,clonal.freq.df[,i+1]))$p.value
  clonal.result.df$odds_ratio[i] <- fisher.test(table(clonal.freq.df$bin,clonal.freq.df[,i+1]))$estimate
}

#compute odds ratio
for(i in 1:nrow(clonal.result.df)){
  tb <-table(clonal.freq.df$bin,clonal.freq.df[,i+1]) + 0.01
  clonal.result.df$odds_ratio[i] <- tb[1,1]*tb[2,2]/(tb[2,1]*tb[1,2])
}

#
clonal.result.df$significance <- ifelse(clonal.result.df$pvalue<0.05,"Yes","No")
clonal.result.df$enrichment <- ifelse(clonal.result.df$odds_ratio < 1,"Concordant low\nenriched","Concordant high\nenriched")
clonal.result.df$enrichment[which(clonal.result.df$significance=="No")] <- "non-significant"

#clonal plot
ggplot(clonal.result.df,aes(log(odds_ratio),-log10(pvalue)))+geom_point(aes(fill=enrichment),pch=21,size=5)+
  scale_x_continuous(expand = c(0,0),limits = c(-2,2),oob = squish)+
  scale_y_continuous(breaks = seq(0,6,1),expand = c(0,0),limits = c(0,6))+
  scale_fill_manual(values = c("#EE0000FF","#3B4992FF","azure4"))+
  ggtitle("Clonal drivers in low vs high")+
  geom_text_repel(data = subset(clonal.result.df,significance=="Yes"),aes(label=Driver),max.overlaps = 25,size=3)+
  labs(fill="")+xlab("Odds ratio")+ylab("-log10(p-value)")+
  geom_vline(xintercept = 0)+geom_hline(yintercept = -log10(0.05),lty="dashed")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))


#subclonal
Subclonal.SNV <- Tx421_regionMut[which(Tx421_regionMut$PyCloneClonal_SC == "S"),c("RegionID","PyCloneClonal_SC","Hugo_Symbol")]
Subclonal.SNV <- as.data.frame(dcast(Subclonal.SNV, RegionID~Hugo_Symbol, value.var = "PyCloneClonal_SC"))
colnames(Subclonal.SNV)[1] <- "sample_name"

#binary df for all primary Tx samples, 0=no mutation, 1=mutation
Subclonal.freq.df <- data.frame(sample_name = gsub("LTX0",Tx421_riskscore$sample_name,replacement = "LTX"), bin = Tx421_riskscore$bin)
Subclonal.freq.df <- left_join(Subclonal.freq.df, Subclonal.SNV)
for(i in 3:ncol(Subclonal.freq.df)){
  Subclonal.freq.df[which(is.na(Subclonal.freq.df[,i])),i] <- 0
}

for(i in 3:ncol(Subclonal.freq.df)){
  Subclonal.freq.df[which(Subclonal.freq.df[,i] > 1),i] <- 1
}

#remove low event number sample
rownames(Subclonal.freq.df) <- Subclonal.freq.df$sample_name
Subclonal.freq.df <- dplyr::select(Subclonal.freq.df, -contains("sample_name"))
Subclonal.freq.df <- Subclonal.freq.df[, -c(which(colSums(Subclonal.freq.df[,-1]) == 0) + 1 )]
Subclonal.freq.df <- Subclonal.freq.df[,-c(which(colSums(Subclonal.freq.df[,-1]) < nrow(Subclonal.freq.df)*0.01)+1)]
Subclonal.freq.df$bin <- factor(Subclonal.freq.df$bin, levels = c("Low","High"))

#Fisher exact test pvalue
Subclonal.result.df <- data.frame(Driver = colnames(Subclonal.freq.df)[-1],pvalue=NA,odds_ratio = NA)

for(i in 1:nrow(Subclonal.result.df)){
  Subclonal.result.df$pvalue[i] <- fisher.test(table(Subclonal.freq.df$bin,Subclonal.freq.df[,i+1]))$p.value
  Subclonal.result.df$odds_ratio[i] <- fisher.test(table(Subclonal.freq.df$bin,Subclonal.freq.df[,i+1]))$estimate
}


#
Subclonal.result.df$significance <- ifelse(Subclonal.result.df$pvalue<0.05,"Yes","No")
Subclonal.result.df$enrichment <- ifelse(Subclonal.result.df$odds_ratio < 1,"Concordant low\nenriched","Concordant high\nenriched")
Subclonal.result.df$enrichment[which(Subclonal.result.df$significance=="No")] <- "non-significant"

#subclonal plot
ggplot(Subclonal.result.df,aes(log(odds_ratio),-log10(pvalue)))+geom_point(aes(fill=enrichment),pch=21,size=5)+
  scale_x_continuous(expand = c(0,0),limits = c(-2,2),oob = squish)+
  scale_y_continuous(breaks = seq(0,6,1),expand = c(0,0),limits = c(0,6))+
  scale_fill_manual(values = c("#EE0000FF","azure4","#3B4992FF"))+
  ggtitle("Subclonal drivers in low vs high")+
  geom_text_repel(data = subset(Subclonal.result.df,significance=="Yes"),aes(label=Driver),max.overlaps = 25,size=3)+
  labs(fill="")+xlab("Odds ratio")+ylab("-log10(p-value)")+
  geom_vline(xintercept = 0)+geom_hline(yintercept = -log10(0.05),lty="dashed")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))


## ED 9B

#label clonal  mutation event
for(i in 1:ncol(clonal.freq.df[,-1])){
  clonal.freq.df[which(clonal.freq.df[,i+1] == 1),i+1] <- "C"
}
for(i in 1:ncol(clonal.freq.df[,-1])){
  clonal.freq.df[which(clonal.freq.df[,i+1] != "C"),i+1] <- NA
}

#add patient-level ORACLE risk-class label
tmp <- clonal.freq.df[,c(1,which(colnames(clonal.freq.df) %in% clonal.result.df$Driver[which(clonal.result.df$pvalue < 0.05)]))]
tmp$sample_name <- gsub("LTX",rownames(tmp),replacement = "LTX0")
tmp <- left_join(tmp, Tx421_riskscore[,c("sample_name", "ORACLE_class")], by = "sample_name")

#order by patient class
tmp$ORACLE_class <- factor(tmp$ORACLE_class, levels = c("Low","Discordant","High"))
rownames(tmp) <- tmp$sample_name
tmp <- dplyr::arrange(tmp, ORACLE_class)

#oncoprint for driver genes enriched for clonal mutation
col = c(C = "#36B084")
color2 <- structure(c("#3B4992FF","azure4","#EE0000FF"), names = c("Low", "Discordant", "High"))
ha <- HeatmapAnnotation("ORACLE" = tmp$ORACLE_class,border = T,show_legend = T,
                        show_annotation_name = F,col = list("ORACLE"=color2))
color3 <- 
  ha2 <- HeatmapAnnotation("ORACLE" = tmp$bin,border = T,show_legend = T,
                           show_annotation_name = F,col = list("ORACLE"=color2))
oncoPrint(t(tmp[,-c(1,12:13)]),
          alter_fun = list( background = alter_graphic("rect", fill = "#DCDCDC"), 
                            C = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.5, gp = gpar(fill = col["C"], col = NA))
                            #mut = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["1"], col = NA))
                            #unknown = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["unknown"], col = NA))
          ),show_column_names = T,show_row_names = T,column_names_gp = gpar(fontsize=4),row_names_gp = gpar(fontsize=7),
          row_names_side = "right",show_pct=T,column_order = tmp$sample_name,top_annotation = ha, 
          col = col)


## 9C

#select significant SCNA in high-risk cohort
HighLow_summary$high_sig <- NA
HighLow_summary$high_sig[which(HighLow_summary$log_qval_high > -log10(0.05) & HighLow_summary$Gscore_diff_raw > 0)] <- HighLow_summary$Type[which(HighLow_summary$log_qval_high > -log10(0.05) & HighLow_summary$Gscore_diff_raw > 0)]
HighLow_summary$high_sig[which(is.na(HighLow_summary$high_sig))] <- "Not significant"

#assign amplification and deletion of cytobands
HighLow_Amp <- HighLow_summary[which(HighLow_summary$Type == "Amp"),]
HighLow_summary$order <- nrow(HighLow_summary):1
HighLow_summary$log_qval_high[which(HighLow_summary$Type == "Del")] <- -HighLow_summary$log_qval_high[which(HighLow_summary$Type == "Del")]

ggplot(HighLow_summary)+
  geom_segment(aes(x=fct_reorder(cytoband_name,order),xend=fct_reorder(cytoband_name,order),y=log_qval_high,yend=0,col=high_sig),alpha=0.9)+
  coord_flip()+xlab("")+
  scale_color_manual(values = c("#EE0000FF","#3B4992FF","azure4"),na.translate=F)+
  geom_vline(xintercept = c(cumsum(rev(table(HighLow_Amp$Chromosome)))+0.5),lty="dashed",col="gray")+
  geom_hline(yintercept = c(log10(0.05),-log10(0.05)),lty="dashed")+
  geom_hline(yintercept = 0)+
  geom_text_repel(data = HighLow_summary[which(HighLow_summary$Gscore_diff_raw > 0 & HighLow_summary$high_sig %in% c("Amp","Del")),],aes(x=fct_reorder(cytoband_name,order),y=log_qval_high,label=cytoband_name),min.segment.length = unit(0.1, "lines"),max.overlaps = 100)+
  scale_y_continuous(name = "q value",expand = c(0,0),breaks = seq(-3,3,1),limits = c(-3.1,3.1))+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.y = element_blank(),axis.ticks.y = element_blank())




### Extended Data Fig 10

#select LUSC in UPPSALA cohort
UPP_LUSC <- UPP_count_data[,c(1,which(colnames(UPP_count_data) %in% colnames(UPP_sample)[which(UPP_sample[11,] == "histology: 1")]))]

#map ensembl ID to gene name
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genesymbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),filters = "ensembl_gene_id",values = UPP_LUSC$Ensembl_gene_id,mart = ensembl)
colnames(UPP_LUSC)[1] <- "ensembl_gene_id"
UPP_LUSC <- left_join(UPP_LUSC, genesymbol, by = "ensembl_gene_id")

#retreive clinical data
UPP_clinical <- UPP_sample[c(10,12,15:16),-1]
rownames(UPP_clinical) <- c("Stage","surgery date","vital date","cens")
UPP_clinical <- as.data.frame( t(UPP_clinical) )

UPP_clinical$`surgery date` <- gsub("surgery date: ", UPP_clinical$`surgery date`, replacement = "")
UPP_clinical$`vital date` <- gsub("vital date: ", UPP_clinical$`vital date`,replacement = "")
UPP_clinical$Stage <- gsub("stage tnm: ", UPP_clinical$Stage, replacement = "")
UPP_clinical$cens <- gsub("dead: ",UPP_clinical$cens,replacement = "")
UPP_clinical <- UPP_clinical[which(rownames(UPP_clinical) %in% colnames(UPP_sample)[which(UPP_sample[11,] == "histology: 1")]),]
UPP_clinical$surv_time <- apply(UPP_clinical, 1, FUN = function(x){ difftime(x[3],x[2]) })
UPP_clinical$cens <- as.numeric(UPP_clinical$cens)

#select max expression for duplicate genes
UPP_LUSC <- UPP_LUSC[which(UPP_LUSC$hgnc_symbol != ""),]
UPP_LUSC <- UPP_LUSC[-which(rowSums(UPP_LUSC[,-c(1, ncol(UPP_LUSC))]) < 10),]
UPP_LUSC <- UPP_LUSC[,-1]
UPP_LUSC <- aggregate(UPP_LUSC[,-ncol(UPP_LUSC)], by = list(Gene = UPP_LUSC$hgnc_symbol), mean)
UPP_LUSC <- column_to_rownames(UPP_LUSC, var = "Gene")

#prepare DESeq Dataset
coldata<-data.frame(sample = colnames(UPP_LUSC),row.names = colnames(UPP_LUSC))
UPP_LUSC <- round(UPP_LUSC)
dds <- DESeqDataSetFromMatrix(countData = UPP_LUSC,
                              colData = coldata,
                              design = ~1)

#vst normalisation
vsd <- vst(dds, blind=TRUE, fitType="mean")
UPP_vsd <- assay(vsd)

#calculate ORACLE RS
UPP_vsd_ORACLE <- UPP_vsd[Supp_Table_5_ORACLE$Gene.Symbol,]
UPP_ORACLE_rs <- data.frame(RiskScore = colSums(UPP_vsd_ORACLE[Supp_Table_5_ORACLE$Gene.Symbol,] * Supp_Table_5_ORACLE$Model.Coefficient), sample = colnames(UPP_vsd_ORACLE))

#join clinical data
UPP_clinical$sample <- rownames(UPP_clinical)
UPP_ORACLE_rs <- left_join(UPP_ORACLE_rs, UPP_clinical, by = "sample")

#censor survival >5 years
UPP_ORACLE_rs$surv_time <- UPP_ORACLE_rs$surv_time/365
UPP_ORACLE_rs$cens[which(UPP_ORACLE_rs$surv_time >5)]<-0
UPP_ORACLE_rs$surv_time[which(UPP_ORACLE_rs$surv_time>5)] <- 5

#UVA cox regression
coxmd <- coxph(Surv(surv_time,cens) ~ RiskScore , data = UPP_ORACLE_rs)
ggforest(coxmd)


