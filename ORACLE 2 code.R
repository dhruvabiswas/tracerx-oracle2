
library(tidyverse)
library(RColorBrewer)
library(survminer)
library(survival)
library(DESeq2)
library(readxl)
library(ggrepel)
library(circlize)
library(maditr)
library(gridExtra)
library(forestplot)
library(fst)
library(stringr)
library(ggalluvial)
library(DescTools)
library(scales)
library(biomaRt)


#load data

all_tumour_df_20220209 <- readRDS("all_tumour_df_20220209.RDS")
all_tumour_df_20220209$tumour_id_mphase_cruk <- apply(all_tumour_df_20220209, 1 , FUN = function(x){ if(grepl("Cluster", x[12])) {x[12] <- gsub(".*-", x[12], replacement = paste0(x[4] , "_"))} else {x[12] <- x[4]} } )

Tx321_riskscore <- read.csv("20230727_TRACERx321_riskscores.csv")

cutoff <- read.csv("2018-11-07_oracle_de.novo_cut.off.csv")
cutoff <- cutoff$RiskScore[1]


### 1B

#remove single region sample
Tx321_samplingbias <- Tx321_riskscore[-which(Tx321_riskscore$High == 0 & Tx321_riskscore$Low == 1),]
Tx321_samplingbias <- Tx321_samplingbias[-which(Tx321_samplingbias$High == 1 & Tx321_samplingbias$Low == 0),]

#join stage info
Tx321_samplingbias <- left_join(Tx321_samplingbias,all_tumour_df_20220209[,c("tumour_id_mphase_cruk","pTNMStage_v8")], by = "tumour_id_mphase_cruk")
Tx321_samplingbias$pTNMStage_v8[which(Tx321_samplingbias$cruk_id == "CRUK0704")] <- "2b" #pathology report - T3N0M0

#group into stage I, II, III
Tx321_samplingbias$pTNMStage_v8 <- gsub("a",Tx321_samplingbias$pTNMStage_v8,replacement = "")
Tx321_samplingbias$pTNMStage_v8 <- gsub("b",Tx321_samplingbias$pTNMStage_v8,replacement = "")

#mean ORACLE RS
for(i in 1:nrow(Tx321_samplingbias)){
  Tx321_samplingbias$ORACLE_mean[i] <- Tx321_samplingbias$ORACLE_riskscore[which(Tx321_samplingbias$cruk_id == Tx321_samplingbias$cruk_id[i])] %>% mean()
}

#arrange order
Tx321_samplingbias$ORACLE_class <- factor(Tx321_samplingbias$ORACLE_class,levels = c("Low","Discordant","High"))
Tx321_samplingbias <- arrange(Tx321_samplingbias,ORACLE_class,ORACLE_mean)
Tx321_samplingbias$order <- 1:nrow(Tx321_samplingbias)

#
pdf("Fig1A_ORACLE_scatter.pdf",width=13,height=3)

ggplot(Tx321_samplingbias,aes(x=fct_reorder(cruk_id, order),y = ORACLE_riskscore))+
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

dev.off()




### 1C

#bin regions into low/high risk
signature_TSB <- Tx321_samplingbias

signature_TSB$Li_bin <- ifelse(signature_TSB$Li_RS < median(signature_TSB$Li_RS), "Low","High")
signature_TSB$Wang_bin <- ifelse(signature_TSB$Wang_RS < median(signature_TSB$Wang_RS), "Low","High")
signature_TSB$Zhao_bin <- ifelse(signature_TSB$Zhao_RS < median(signature_TSB$Zhao_RS), "Low","High")
signature_TSB$Song_bin <- ifelse(signature_TSB$Song_RS < median(signature_TSB$Song_RS), "Low","High")
signature_TSB$Jin_bin <- ifelse(signature_TSB$Jin_RS < median(signature_TSB$Jin_RS), "Low","High")
signature_TSB$LiFeng_bin <- ifelse(signature_TSB$LiFeng_RS < median(signature_TSB$LiFeng_RS), "Low","High")

#
tmp <- table(signature_TSB$cruk_id,signature_TSB$Li_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Li_class <- NA
tmp$Li_class <- ifelse(tmp$High > 0, paste(tmp$Li_class, "High", sep=""), tmp$Li_class)
tmp$Li_class <- ifelse(tmp$Low > 0, paste(tmp$Li_class, "Low", sep=""), tmp$Li_class)
tmp$Li_class <- gsub(x=tmp$Li_class, pattern="NA", replacement="")
tmp$Li_class <- gsub(x=tmp$Li_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Li_class <- factor(tmp$Li_class, levels = c("Low","Discordant","High"))
signature_TSB <- left_join(signature_TSB, tmp[,c("cruk_id","Li_class")], by = "cruk_id")

#
tmp <- table(signature_TSB$cruk_id,signature_TSB$Wang_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Wang_class <- NA
tmp$Wang_class <- ifelse(tmp$High > 0, paste(tmp$Wang_class, "High", sep=""), tmp$Wang_class)
tmp$Wang_class <- ifelse(tmp$Low > 0, paste(tmp$Wang_class, "Low", sep=""), tmp$Wang_class)
tmp$Wang_class <- gsub(x=tmp$Wang_class, pattern="NA", replacement="")
tmp$Wang_class <- gsub(x=tmp$Wang_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Wang_class <- factor(tmp$Wang_class, levels = c("Low","Discordant","High"))
signature_TSB <- left_join(signature_TSB, tmp[,c("cruk_id","Wang_class")], by = "cruk_id")

#
tmp <- table(signature_TSB$cruk_id,signature_TSB$Zhao_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Zhao_class <- NA
tmp$Zhao_class <- ifelse(tmp$High > 0, paste(tmp$Zhao_class, "High", sep=""), tmp$Zhao_class)
tmp$Zhao_class <- ifelse(tmp$Low > 0, paste(tmp$Zhao_class, "Low", sep=""), tmp$Zhao_class)
tmp$Zhao_class <- gsub(x=tmp$Zhao_class, pattern="NA", replacement="")
tmp$Zhao_class <- gsub(x=tmp$Zhao_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Zhao_class <- factor(tmp$Zhao_class, levels = c("Low","Discordant","High"))
signature_TSB <- left_join(signature_TSB, tmp[,c("cruk_id","Zhao_class")], by = "cruk_id")

#
tmp <- table(signature_TSB$cruk_id,signature_TSB$Song_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Song_class <- NA
tmp$Song_class <- ifelse(tmp$High > 0, paste(tmp$Song_class, "High", sep=""), tmp$Song_class)
tmp$Song_class <- ifelse(tmp$Low > 0, paste(tmp$Song_class, "Low", sep=""), tmp$Song_class)
tmp$Song_class <- gsub(x=tmp$Song_class, pattern="NA", replacement="")
tmp$Song_class <- gsub(x=tmp$Song_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Song_class <- factor(tmp$Song_class, levels = c("Low","Discordant","High"))
signature_TSB <- left_join(signature_TSB, tmp[,c("cruk_id","Song_class")], by = "cruk_id")

#
tmp <- table(signature_TSB$cruk_id,signature_TSB$Jin_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$Jin_class <- NA
tmp$Jin_class <- ifelse(tmp$High > 0, paste(tmp$Jin_class, "High", sep=""), tmp$Jin_class)
tmp$Jin_class <- ifelse(tmp$Low > 0, paste(tmp$Jin_class, "Low", sep=""), tmp$Jin_class)
tmp$Jin_class <- gsub(x=tmp$Jin_class, pattern="NA", replacement="")
tmp$Jin_class <- gsub(x=tmp$Jin_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$Jin_class <- factor(tmp$Jin_class, levels = c("Low","Discordant","High"))
signature_TSB <- left_join(signature_TSB, tmp[,c("cruk_id","Jin_class")], by = "cruk_id")

#
tmp <- table(signature_TSB$cruk_id,signature_TSB$LiFeng_bin)
tmp <- data.frame(High=as.matrix(tmp)[,"High"], Low=as.matrix(tmp)[,"Low"])
tmp$LiFeng_class <- NA
tmp$LiFeng_class <- ifelse(tmp$High > 0, paste(tmp$LiFeng_class, "High", sep=""), tmp$LiFeng_class)
tmp$LiFeng_class <- ifelse(tmp$Low > 0, paste(tmp$LiFeng_class, "Low", sep=""), tmp$LiFeng_class)
tmp$LiFeng_class <- gsub(x=tmp$LiFeng_class, pattern="NA", replacement="")
tmp$LiFeng_class <- gsub(x=tmp$LiFeng_class, pattern="HighLow", replacement="Discordant")
tmp$cruk_id <- rownames(tmp)
tmp$LiFeng_class <- factor(tmp$LiFeng_class, levels = c("Low","Discordant","High"))
signature_TSB <- left_join(signature_TSB, tmp[,c("cruk_id","LiFeng_class")], by = "cruk_id")

write.csv(signature_TSB, "signature_TSB.csv", row.names = F)

signature_TSB_summary <- select(signature_TSB, cruk_id, ORACLE_class, Li_class, Wang_class, Zhao_class, Song_class, Jin_class, LiFeng_class)
signature_TSB_summary <- signature_TSB_summary[which(!duplicated(signature_TSB_summary$cruk_id)),]


#
TSB <- rbind(table(signature_TSB_summary$ORACLE_class)*100/nrow(signature_TSB_summary),
             table(signature_TSB_summary$Li_class)*100/nrow(signature_TSB_summary),
             table(signature_TSB_summary$Wang_class)*100/nrow(signature_TSB_summary),
             table(signature_TSB_summary$Zhao_class)*100/nrow(signature_TSB_summary),
             table(signature_TSB_summary$Song_class)*100/nrow(signature_TSB_summary),
             table(signature_TSB_summary$Jin_class)*100/nrow(signature_TSB_summary),
             table(signature_TSB_summary$LiFeng_class)*100/nrow(signature_TSB_summary)) %>% as.data.frame()

#tidy up variable level
TSB$signature <- c("ORACLE","Li JAMA Oncol 2017","Wang Front Immunol 2022","Zhao Lung Cancer 2020","Song Sci Rep 2022","Jin J Immunol Res 2022","Li Sci Rep 2022")
TSB$signature <- factor(TSB$signature, levels = c("ORACLE" ,"Li JAMA Oncol 2017","Wang Front Immunol 2022","Zhao Lung Cancer 2020","Song Sci Rep 2022","Jin J Immunol Res 2022","Li Sci Rep 2022"))

TSB_melt <- melt(TSB)
TSB_melt$variable <- factor(TSB_melt$variable, levels = c("Low","High","Discordant"))

pdf("Fig1B_TSB_pie.pdf",width=5,height=5)
ggplot(TSB_melt[which(TSB_melt$signature == "ORACLE"),],aes(x = signature,y=value,fill=variable))+
  geom_col(width = 0.6,position = "stack", col = "black")+
  coord_polar("y", start=0)+
  ylab("")+ xlab("")+ ggtitle("ORACLE") +
  geom_text(aes( x=signature, y = value-0.5, label = paste0(signif(value,digits = 2),"%")),position=position_stack(vjust=0.5),size=8)+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme_void()+
  theme(legend.position = "none",)

ggplot(TSB_melt[which(TSB_melt$signature == "Li JAMA Oncol 2017"),],aes(x = signature,y=value,fill=variable))+
  geom_col(width = 0.6,position = "stack", col = "black")+
  coord_polar("y", start=0)+
  ylab("")+ xlab("")+ ggtitle("Li JAMA Oncol 2017") +
  geom_text(aes( x=signature, label = paste0(signif(value,digits = 2),"%")),position=position_stack(vjust=0.5),size=8)+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme_void()+
  theme(legend.position = "none",)

ggplot(TSB_melt[which(TSB_melt$signature == "Wang Front Immunol 2022"),],aes(x = signature,y=value,fill=variable))+
  geom_col(width = 0.6,position = "stack", col = "black")+
  coord_polar("y", start=0)+
  ylab("")+ xlab("")+ ggtitle("Wang Front Immunol 2022") +
  geom_text(aes( x=signature, y = value, label = paste0(signif(value,digits = 2),"%")),position=position_stack(vjust=0.5),size=8)+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme_void()+
  theme(legend.position = "none",)

ggplot(TSB_melt[which(TSB_melt$signature == "Zhao Lung Cancer 2020"),],aes(x = signature,y=value,fill=variable))+
  geom_col(width = 0.6,position = "stack", col = "black")+
  coord_polar("y", start=0)+
  ylab("")+ xlab("")+ ggtitle("Zhao Lung Cancer 2020") +
  geom_text(aes( x=signature, y = value, label = paste0(signif(value,digits = 2),"%")),position=position_stack(vjust=0.5),size=8)+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme_void()+
  theme(legend.position = "none",)

ggplot(TSB_melt[which(TSB_melt$signature == "Song Sci Rep 2022"),],aes(x = signature,y=value,fill=variable))+
  geom_col(width = 0.6,position = "stack", col = "black")+
  coord_polar("y", start=0)+
  ylab("")+ xlab("")+ ggtitle("Song Sci Rep 2022") +
  geom_text(aes( x=signature, y = value, label = paste0(signif(value,digits = 2),"%")),position=position_stack(vjust=0.5),size=8)+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme_void()+
  theme(legend.position = "none",)

ggplot(TSB_melt[which(TSB_melt$signature == "Jin J Immunol Res 2022"),],aes(x = signature,y=value,fill=variable))+
  geom_col(width = 0.6,position = "stack", col = "black")+
  coord_polar("y", start=0)+
  ylab("")+ xlab("")+ ggtitle("Jin J Immunol Res 2022") +
  geom_text(aes( x=signature, y = value, label = paste0(signif(value,digits = 2),"%")),position=position_stack(vjust=0.5),size=8)+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme_void()+
  theme(legend.position = "none",)

ggplot(TSB_melt[which(TSB_melt$signature == "Li Sci Rep 2022"),],aes(x = signature,y=value,fill=variable))+
  geom_col(width = 0.6,position = "stack", col = "black")+
  coord_polar("y", start=0)+
  ylab("")+ xlab("")+ ggtitle("Li Sci Rep 2022") +
  geom_text(aes( x=signature, y = value, label = paste0(signif(value,digits = 2),"%")),position=position_stack(vjust=0.5),size=8)+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF","azure4"))+
  theme_void()+
  theme(legend.position = "none",)

dev.off()

write.csv(TSB, "Discordant_perc.csv", row.names = F)




### 1D

#load data
TSB_metric1 <- read.csv("Discordant_perc.csv")
TSB_metric2 <- read.csv("Cluster_concordance_AUC.csv")
TSB_metric3 <- read.csv("ExpVar_byHouseham.csv")
TSB_metric4 <- read.csv("LeastBiopsy_TSB.csv")

#tidy up
TSB_metric1 <- select(TSB_metric1, signature, Discordant)
TSB_metric2 <- select(TSB_metric2, signature, AUC)
TSB_metric3 <- aggregate(TSB_metric3$mean_var, by = list(signature = TSB_metric3$signature), median)

#ranking by each metric
TSB_metric1 <- arrange(TSB_metric1, Discordant)
TSB_metric1$rank <- 1:nrow(TSB_metric1)

TSB_metric2 <- arrange(TSB_metric2, desc(AUC))
TSB_metric2$rank <- 1:nrow(TSB_metric2)

TSB_metric3 <- arrange(TSB_metric3, x)
TSB_metric3$rank <- 1:nrow(TSB_metric3)

TSB_metric4 <- arrange(TSB_metric4, biop_need)
TSB_metric4$rank <- 1:nrow(TSB_metric4)

#
colnames(TSB_metric1)[2] <- "value"
colnames(TSB_metric2)[2] <- "value"
colnames(TSB_metric3)[2] <- "value"
colnames(TSB_metric4)[2] <- "value"

#
TSB_metric1$value <- paste0(signif(TSB_metric1$value,digits = 2), "%")
TSB_metric2$value <- signif(TSB_metric2$value,digits = 3)
TSB_metric3$value <- signif(TSB_metric3$value,digits = 2)
TSB_metric4$value <- signif(TSB_metric4$value,digits = 3)

#
TSB_metric1 <- select(TSB_metric1,signature, rank,value)
TSB_metric2 <- select(TSB_metric2,signature, rank,value)
TSB_metric3 <- select(TSB_metric3,signature, rank,value)
TSB_metric4 <- select(TSB_metric4,signature, rank,value)

#assign metric names
TSB_metric1$metric <- "Discordant rate"
TSB_metric2$metric <- "AUC"
TSB_metric3$metric <- "Expression SD"
TSB_metric4$metric <- "Least Biopsy"

All_TSB_metrics <- rbind(TSB_metric1, TSB_metric2, TSB_metric3, TSB_metric4)

#
tmp <- aggregate(All_TSB_metrics$rank , by = list(signature = All_TSB_metrics$signature), mean)
tmp <- arrange(tmp, x)
tmp$rank <- 1:nrow(tmp)
colnames(tmp)[2] <- "value"
tmp$metric <- "mean"
tmp <- select(tmp, signature, rank, value, metric)
tmp$rank[which(tmp$value == 4.5)] <- 3

#
All_TSB_metrics <- rbind(All_TSB_metrics, tmp)
All_TSB_metrics$rank <- factor(All_TSB_metrics$rank,levels = c(7,6,5,4,3,2,1))
All_TSB_metrics$metric <- factor(All_TSB_metrics$metric,levels = c("Discordant rate","AUC","Expression SD","Least Biopsy","mean"))
All_TSB_metrics$signature <- factor(All_TSB_metrics$signature, levels = c("ORACLE","Li JAMA Oncol 2017","Wang Front Immunol 2022","Zhao Lung Cancer 2020","Song Sci Rep 2022","Jin J Immunol Res 2022","Li Sci Rep 2022"))
All_TSB_metrics$label <- ifelse(All_TSB_metrics$metric == "mean", "Yes","No")


pdf("All_TSBmetrics_rank.pdf",width = 12,height=3)
ggplot(All_TSB_metrics, aes(metric, rank)) + 
  geom_point(aes(col=signature,size = label))+
  geom_line(aes(group = signature,col=signature),size=2) +
  xlab("") + ylab("Ranking") +
  geom_text(aes(x=metric, y=rank, label = value)) + 
  scale_size_manual(values = c(5,10)) +
  scale_color_manual(values = c("#A6CEE3","#1F78B4","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")) +
  theme(panel.background = element_blank(),panel.border = element_blank(),axis.ticks = element_blank())
dev.off()

write.csv(All_TSB_metrics, "TSB_fourMetrics_results.csv", row.names = F)






#load data

Tx321_riskscore <- read.csv("20230727_TRACERx321_riskscores.csv")
Tx421_riskscore <- read.csv("20230727_TRACERx421_riskscores.csv")

TRACERx_clinical <- readRDS("all_patient_df_20221109.RDS")

all_tumour_df_20220209 <- readRDS("all_tumour_df_20220209.RDS")
all_tumour_df_20220209$tumour_id_mphase_cruk <- apply(all_tumour_df_20220209, 1 , FUN = function(x){ if(grepl("Cluster", x[12])) {x[12] <- gsub(".*-", x[12], replacement = paste0(x[4] , "_"))} else {x[12] <- x[4]} } )





### 2A

#overall survival data
Tx321_Survdata <- left_join(Tx321_riskscore, TRACERx_clinical[,c("cruk_id","cens_os","os_time","lung_specific_time","cens_lung_specific")], by = "cruk_id")

#join clinicopathological
Tx321_Survdata <- left_join(Tx321_Survdata, all_tumour_df_20220209[,c("tumour_id_mphase_cruk","age","sex","pack_years_calculated","pTNMStage_v8","adjuvant_treatment_YN")], by = "tumour_id_mphase_cruk")
Tx321_Survdata$pTNMStage_v8[which(grepl("CRUK0704",Tx321_Survdata$tumour_id_mphase_cruk))] <- "2b" #collision tumour

#censor survival >5 years
#OS
Tx321_Survdata$os_time <- Tx321_Survdata$os_time/365
Tx321_Survdata$cens_os[which(Tx321_Survdata$os_time >5)]<-0
Tx321_Survdata$os_time[which(Tx321_Survdata$os_time>5)] <- 5

#LSS
Tx321_Survdata$lung_specific_time <- Tx321_Survdata$lung_specific_time/365
Tx321_Survdata$cens_lung_specific[which(Tx321_Survdata$lung_specific_time >5)]<-0
Tx321_Survdata$lung_specific_time[which(Tx321_Survdata$lung_specific_time>5)] <- 5

#DFS

#
Tx321_tumour_Survdata <- Tx321_Survdata[which(!duplicated(Tx321_Survdata$cruk_id)),]
Tx321_tumour_Survdata$ORACLE_class <- factor(Tx321_tumour_Survdata$ORACLE_class, levels = c("Low","Discordant","High"))

#KM plot
j <- ggsurvplot(fit = survfit(Surv(os_time,cens_os)~ORACLE_class,data = Tx321_tumour_Survdata),
                legend.labs = c( "Low","Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                title="TRACERx validation cohort",
                legend="none", xlab="Time (years)", ylab="Overall Survival",
                pval = T,pval.size = 6,pval.coord = c(0, 0.05),
                ggtheme = theme(title = element_text(size = 18),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j

#hazard ratio
coxmd <- coxph(Surv(os_time, cens_os) ~ ORACLE_class , data = Tx321_tumour_Survdata)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])





### 2B

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

#MVA
MVA.cox.res <- coxph(formula = Surv(Tx321_tumour_Survdata$os_time, Tx321_tumour_Survdata$cens_os) ~ 
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

#
MVA.cox.res$pvalue <-signif(MVA.cox.res$P,digits = 1) 
MVA.cox.res$pvalue[which(is.na(MVA.cox.res$pvalue))] <- "Ref"
MVA.cox.res$mean <- round(MVA.cox.res$HR,digits = 2)

for(i in 1:nrow(MVA.cox.res)){
  lci <- round(MVA.cox.res$lower_ci[i],digits=2);
  uci <- round(MVA.cox.res$upper_ci[i],digits=2);
  MVA.cox.res$result[i] <- paste(lci,uci,sep = " - ")
}

#
MVA.cox.res$mean[which(MVA.cox.res$pvalue == "Ref")] <- ""

#
MVA.cox.res <- rbind(c(1,NA,NA,"","","p-value","HR","CI"),MVA.cox.res)

#turn data from chr to numeric
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

###Forest plot
pdf("Forestplot_RNAseq_sig.pdf",width=11,height=5)
forestplot(MVA.cox.res[,c(5:8)],graph.pos=2,
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





### 2C

sample.df  <- Tx321_riskscore
colnames(sample.df)[5] <- "ORACLE_RS"


#simulate sampling 1 region per tumour
create_surv <- function(signature){
  
  #
  RiskScore.df <- select(sample.df, cruk_id, paste(signature,"RS",sep = "_"))
  
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


#survival association - log-rank p, cox regression
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


ORACLE.cox.df <- Reduce(rbind,ORACLE.cox.df)
Li.cox.df <- Reduce(rbind,Li.cox.df)
Wang.cox.df <- Reduce(rbind,Wang.cox.df)
Zhao.cox.df <- Reduce(rbind,Zhao.cox.df)
Song.cox.df <- Reduce(rbind,Song.cox.df)
Jin.cox.df <- Reduce(rbind,Jin.cox.df)
LiFeng.cox.df <- Reduce(rbind,LiFeng.cox.df)

##boot HR
#
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

#
ggplot(data=ORACLE.cox.df,aes(HR))+
  geom_ribbon(aes(xmax = mean(cox_lci), xmin = mean(cox_uci),ymax = 3, ymin= 0),alpha=0.2, fill = "#FDAE5F") + 
  geom_density(fill="#4EB0EC",alpha=0.4)+ 
  xlab("Hazard ratio") + ylab("Density")+
  scale_x_continuous(breaks = seq(0,3.5,1),expand = c(0,0),limits = c(0,3.5))+
  scale_y_continuous(breaks = seq(0,3,1),expand = c(0,0),limits = c(0,3))+
  geom_vline(xintercept = 1, lty="dashed",col="#25508F")+
  geom_vline(xintercept = ORACLE_HR_mu, col = "red") +
  theme(plot.margin = unit(c(2,4,2,2),"mm"),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))





### 2D 
#select stage I patients
Survdata_stageI <- Tx321_tumour_Survdata[which(Tx321_tumour_Survdata$pTNMStage_v8 == "1a" | Tx321_tumour_Survdata$pTNMStage_v8 == "1b"),]

#KM plot by clinical stage
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

#Hazard ratio
coxmd <- coxph(Surv(os_time, cens_os) ~ pTNMStage_v8 , data = Survdata_stageI)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])


#KM plot classified by ORACLE
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

#hazard ratio
coxmd <- coxph(Surv(os_time, cens_os) ~ ORACLE_class , data = Survdata_stageI)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])





### 3A
#KM plot - lung cancer specific survival
j <- ggsurvplot(fit = survfit(Surv(lung_specific_time,cens_lung_specific)~ORACLE_class,data = Tx321_tumour_Survdata),
                legend.labs = c( "Low","Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                title="TRACERx validation cohort",
                legend="none", xlab="Time (years)", ylab="Lung-specific Survival",
                pval = T,pval.size = 6,pval.coord = c(0, 0.05),
                ggtheme = theme(title = element_text(size = 18),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j

#hazard ratio
coxmd <- coxph(Surv(lung_specific_time, cens_lung_specific) ~ ORACLE_class , data = Tx321_tumour_Survdata)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])





### 4C

#split into adjuvant/no adjuvant cohorts
Survdata_No_adjuvant <- Tx321_tumour_Survdata[which( Tx321_tumour_Survdata$adjuvant_treatment_YN == "No adjuvant"),]
Survdata_with_adjuvant <- Tx321_tumour_Survdata[which( Tx321_tumour_Survdata$adjuvant_treatment_YN == "Adjuvant"),]

#KM plot
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

#hazard ratio
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

#hazard ratio
coxmd <- coxph(Surv(os_time, cens_os) ~ ORACLE_class , data = Survdata_with_adjuvant)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])





### ED 5A
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





### ED 5B
#MVA
MVA.cox.res <- coxph(formula = Surv(Tx321_tumour_Survdata$lung_specific_time, Tx321_tumour_Survdata$cens_lung_specific) ~ 
                       Tx321_tumour_Survdata$ORACLE_mean + 
                       Tx321_tumour_Survdata$sex + 
                       Tx321_tumour_Survdata$age + 
                       Tx321_tumour_Survdata$pack_years_calculated + 
                       Tx321_tumour_Survdata$adjuvant_treatment_YN +
                       Tx321_tumour_Survdata$TNM_combine )

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

#
MVA.cox.res$pvalue <-signif(MVA.cox.res$P,digits = 1) 
MVA.cox.res$pvalue[which(is.na(MVA.cox.res$pvalue))] <- "Ref"
MVA.cox.res$mean <- round(MVA.cox.res$HR,digits = 2)

for(i in 1:nrow(MVA.cox.res)){
  lci <- round(MVA.cox.res$lower_ci[i],digits=2);
  uci <- round(MVA.cox.res$upper_ci[i],digits=2);
  MVA.cox.res$result[i] <- paste(lci,uci,sep = " - ")
}

#
MVA.cox.res$mean[which(MVA.cox.res$pvalue == "Ref")] <- ""

#
MVA.cox.res <- rbind(c(1,NA,NA,"","","p-value","HR","CI"),MVA.cox.res)

#turn data from chr to numeric
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#Forest plot - lung cancer specific survival
forestplot(MVA.cox.res[,c(5:8)],graph.pos=3,
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






### ED 5C
#lung cancer specific survival KM plot by clinical stage
j <- ggsurvplot(fit = survfit(Surv(lung_specific_time,cens_lung_specific)~pTNMStage_v8,data = Survdata_stageI),
                legend.labs = c("IA", "IB"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","#EE0000FF"),size = 1.5, 
                title="Substaging - StageI", 
                legend="none", xlab="Time (Years)", ylab="Lung-specific survival",
                pval = T,pval.size = 8,
                ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j

#Hazard ratio
coxmd <- coxph(Surv(lung_specific_time,cens_lung_specific) ~ pTNMStage_v8 , data = Survdata_stageI)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])


#KM plot classified by ORACLE
Survdata_stageI$ORACLE_class <- factor(Survdata_stageI$ORACLE_class,levels = c("Low","Discordant","High"))
j<-ggsurvplot(fit = survfit(Surv(lung_specific_time,cens_lung_specific)~ORACLE_class,data = Survdata_stageI),
              legend.labs = c("Low", "Discordant","High"),
              risk.table = TRUE, 
              palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
              title="ORACLE", 
              legend="none", xlab="Time (Years)", ylab="Lung-specific Survival",
              pval = T,pval.size = 8,
              ggtheme = theme(title = element_text(size = 20),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                              panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
              tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j

#hazard ratio
coxmd <- coxph(Surv(lung_specific_time,cens_lung_specific) ~ ORACLE_class , data = Survdata_stageI)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])






### ED 7A

#split by node status
Survdata_node_neg <- Tx321_tumour_Survdata[which(Tx321_tumour_Survdata$pN_stage == "0" | Tx321_tumour_Survdata$pN_stage == "X"),]
Survdata_node_pos <- Tx321_tumour_Survdata[which(Tx321_tumour_Survdata$pN_stage != "0" & Tx321_tumour_Survdata$pN_stage != "X"),]

#KM plot
#Node negative without AD
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

#Node negative with AD
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


#Node positive without AD
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


#Node positive with AD
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






##load data

#Mascaux dataset
Mascaux <- read.delim("GSE33479.txt.gz", as.is=T, check.names=FALSE)
gene.annot <- read.delim("GPL6480-9577.txt", as.is=T, check.names=FALSE,skip=17)
patient.annot <- read.table("GSE33479_series_matrix.txt.gz", header = TRUE,fill = TRUE,skip=25)

Supp_Table_5_ORACLE <- read.csv("Supp_Table_5_ORACLE.csv")

#met clone dataset
Tx321_riskscore <- read.csv("20230727_TRACERx321_riskscores.csv")
Tx421_riskscore <- read.csv("20230727_TRACERx421_riskscores.csv")
Tx421_full_riskscore <- read.csv("20230727_TRACERx421_wholecohort_riskscores.csv")

TRACERx_clinical <- readRDS("all_patient_df_20221109.RDS")
Tx421_clinicopatho <- read.csv("Tx421_merged_clinicopathological_data_20220726.nooutcome.csv")
all_tumour_df_20220209 <- readRDS("all_tumour_df_20220209.RDS")
all_tumour_df_20220209$tumour_id_mphase_cruk <- apply(all_tumour_df_20220209, 1 , FUN = function(x){ if(grepl("Cluster", x[12])) {x[12] <- gsub(".*-", x[12], replacement = paste0(x[4] , "_"))} else {x[12] <- x[4]} } )

seed.region <- read.table("seedingRegionInfo.txt",header = T,sep = "\t")

#CCLE dataset
IC50_results <- read.csv("CCLE_result_table_20230308.csv",row.names = 1)




### 3B

#remove non-matched probes
Mascaux <- Mascaux[-which(is.na(Mascaux$ID_REF)),]
rownames(Mascaux) <- Mascaux$ID_REF
Mascaux <- Mascaux[,-1]

#match gene names to probes
Mascaux$ID <- rownames(Mascaux)
Mascaux <- left_join(Mascaux,gene.annot[,c("ID","GENE_SYMBOL")],by="ID")
Mascaux <- Mascaux[-which(Mascaux$GENE_SYMBOL == ""),]
Mascaux <- Mascaux[,-c(123:127)]

#paste patient ID
colnames(patient.annot) <- paste(patient.annot[15,], colnames(patient.annot),sep = "-")

#duplicated genes - median value
Mascaux_ORACLE <- Mascaux[which(Mascaux$GENE_SYMBOL %in% Supp_Table_5_ORACLE$Gene.Symbol),]
tmp <- aggregate(Mascaux_ORACLE[,-which(colnames(Mascaux_ORACLE) == "GENE_SYMBOL")],by=list(Gene = Mascaux_ORACLE$GENE_SYMBOL),max)

Mascaux_ORACLE <- tmp
rownames(Mascaux_ORACLE) <- Mascaux_ORACLE$Gene
Mascaux_ORACLE <- Mascaux_ORACLE[,-1]

#
Mascaux_ORACLE_RS <- data.frame(riskscore = colSums(Mascaux_ORACLE[Supp_Table_5_ORACLE$Gene.Symbol,]*Supp_Table_5_ORACLE$Model.Coefficient))

#
rownames(Mascaux_ORACLE_RS) <- colnames(patient.annot)[-1]

#
Mascaux_ORACLE_RS$patient <- gsub("-.*", rownames(Mascaux_ORACLE_RS),replacement = "")

#group
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

#compute linear mixed-effect model to compare between pairs of stages
Mascaux_ORACLE_RS$group <- factor(Mascaux_ORACLE_RS$group, levels = c("normal","hyperplasia","metaplasia","mild\ndysplasia","moderate\ndysplasia","severe\ndysplasia","CIS","SCC"))
md1 <- nlme::lme(data=Mascaux_ORACLE_RS[which(Mascaux_ORACLE_RS$group %in% c("normal","metaplasia")),], riskscore ~ group, random = ~1|patient)
r1 <- summary(md1)

md2 <- nlme::lme(data=Mascaux_ORACLE_RS[which(Mascaux_ORACLE_RS$group %in% c("metaplasia","SCC")),], riskscore ~ group, random = ~1|patient)
r2 <- summary(md2)

md3 <- nlme::lme(data=Mascaux_ORACLE_RS[which(Mascaux_ORACLE_RS$group %in% c("metaplasia","hyperplasia")),], riskscore ~ group, random = ~1|patient)
r3 <- summary(md3)

r <- data.frame(group = c("metaplasia","SCC"), pval = c(r1$tTable[2,5],r2$tTable[2,5]))

ggplot(Mascaux_ORACLE_RS,aes(group,riskscore))+geom_boxplot(aes(fill=group))+
  xlab("Developmental stages")+ylab("ORACLE Risk score")+
  scale_y_continuous(expand = c(0,0),limits = c(0,3))+
  scale_fill_manual(values = brewer.pal(8,"Oranges"))+
  geom_text(data = r ,aes(x = group, y = 2.8, label = signif(pval,digits = 2))) +
  theme(legend.position="none",panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))





### 3C

Tx421_full_riskscore <- left_join(Tx421_full_riskscore, Tx421_clinicopatho[,c("sample_name","sample_name_cruk","tumour_id_mphase_cruk","cruk_id")], by = "sample_name")

Rec.sample.df <- Tx421_full_riskscore[which(Tx421_full_riskscore$tumour_id_mphase_cruk %in% Tx421_riskscore$tumour_id_mphase_cruk),]
Rec.sample.df <- Rec.sample.df[which(grepl("SU_T", Rec.sample.df$sample_name )),] #select primary sample
Rec.sample.df$PrimaryRegion <- gsub("-", Rec.sample.df$sample_name_cruk, replacement = ".")

#add metastatic seeding/non-seeding info
Rec.sample.df <- left_join(Rec.sample.df, seed.region[,c("PrimaryRegion","Metastasizing")], by ="PrimaryRegion")
Rec.sample.df <- Rec.sample.df[which(!is.na(Rec.sample.df$Metastasizing)),]

#keep samples with at least one seeding and one non-seeding
remove <- table(Rec.sample.df$cruk_id, Rec.sample.df$Metastasizing) %>% as.data.frame()
Rec.sample.df <- Rec.sample.df[-which(Rec.sample.df$cruk_id %in% remove$Var1[which(remove$Freq == 0)]),]

#select recurrence samples
Rec.sample.df <- left_join(Rec.sample.df, TRACERx_clinical[,c("cruk_id","Relapse_cat_new","Relapse_cat")], by = "cruk_id")
Rec.sample.df <- Rec.sample.df[which(Rec.sample.df$Relapse_cat_new != "No rec"),]

#
Rec.sample.df$Metastasizing[which(Rec.sample.df$Metastasizing == "TRUE")] <- "Seeding"
Rec.sample.df$Metastasizing[which(Rec.sample.df$Metastasizing == "FALSE")] <- "Non-seeding"

#compute linear mixed-effect model, comparing ORACLE score between seeding groups
md <-nlme::lme(data=Rec.sample.df,ORACLE_riskscore~Metastasizing,na.action = na.omit,random=~1|cruk_id)
r <- summary(md)

ggplot(Rec.sample.df,aes(Metastasizing,ORACLE_riskscore))+
  geom_boxplot(aes(fill=Metastasizing), width=0.5)+
  geom_jitter(width = 0.2)+
  ylab("ORACLE Risk-score")+xlab("")+
  scale_fill_manual(values=c("#7DB5FF","#FDCF0B"))+
  scale_y_continuous(breaks=seq(9,12,1),expand = c(0,0),limits = c(9,12))+
  geom_text(aes(x=2,y=11.8,label=paste0("p = ",signif(r$tTable[2,5],digits = 2))),size=4) +
  theme(legend.position = "none",panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())






### 3D

#survival data for validation cohort
Relapse_survdata <- dplyr::select(Tx321_riskscore, ORACLE_class, ORACLE_riskscore,cruk_id,tumour_id_mphase_cruk)
tmp <- aggregate(Relapse_survdata$ORACLE_riskscore, by = list(cruk_id = Relapse_survdata$cruk_id), mean)
colnames(tmp)[2] <- "ORACLE_mean"
Relapse_survdata <- left_join(Relapse_survdata, tmp, by = "cruk_id")
Relapse_survdata <- Relapse_survdata[which(!duplicated(Relapse_survdata$cruk_id)),]

#days to years
Relapse_survdata <- left_join(Relapse_survdata, TRACERx_clinical[,c("cruk_id","dfs_time","cens_dfs","dfs_time_any_event","Relapse_cat_new","first_dfs_any_event")], by = "cruk_id")
Relapse_survdata$dfs_time_any_event <- Relapse_survdata$dfs_time_any_event/365
Relapse_survdata$dfs_time <- Relapse_survdata$dfs_time/365

# For DFS, censor 4 patients who are currently marked as "recurrence" but uncertain whether the recurrence is from 1st primary or 2nd primary (if from 2nd primary, these cases should NOT be marked as recurrence in TRACERx protocol)
#  "CRUK0512", "CRUK0373","CRUK0428","CRUK0511" : currently marked as recurrence, after 2nd primary cancer was confirmed
Relapse_survdata$cens_dfs[which(Relapse_survdata$cruk_id %in% c("CRUK0512","CRUK0373","CRUK0428","CRUK0511"))] <- 0
Relapse_survdata$dfs_time[which(Relapse_survdata$cruk_id %in% c("CRUK0512","CRUK0373","CRUK0428","CRUK0511"))] <- Relapse_survdata$dfs_time_any_event[which(Relapse_survdata$cruk_id %in% c("CRUK0512","CRUK0373","CRUK0428","CRUK0511"))]

#censor patients with survival >5 years
Relapse_survdata$cens_dfs[which(Relapse_survdata$dfs_time >5)]<-0
Relapse_survdata$dfs_time[which(Relapse_survdata$dfs_time >5)] <- 5

#disease-free survival to predict relapse time
Relapse_survdata$ORACLE_class <- factor(Relapse_survdata$ORACLE_class , levels = c("Low","Discordant","High"))
j <- ggsurvplot(fit = survfit(Surv(dfs_time,cens_dfs)~ORACLE_class,data = Relapse_survdata),
                legend.labs = c( "Low","Discordant","High"),
                risk.table = TRUE, 
                palette = c("#3B4992FF","azure4","#EE0000FF"),size = 1.5, 
                title="",
                #ncensor.plot = TRUE,
                legend="none", xlab="Time (years)", ylab="Disease-Free Survival",
                pval = T,pval.size = 6,pval.coord = c(0, 0.05),
                ggtheme = theme(title = element_text(size = 18),axis.text = element_text(size = 18),axis.title = element_text(size = 20),
                                panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank()),
                tables.theme = theme_cleantable())

j$plot <- j$plot+scale_y_continuous(breaks = seq(0,1,0.2))
j$table$theme$panel.border <- element_blank()
j$table$theme$title <- element_text(size = 16)
j

#hazard ratio
coxmd <- coxph(Surv(dfs_time, cens_dfs) ~ ORACLE_class , data = Relapse_survdata)
cox_result <- summary(coxmd)
cox_result <- data.frame(HR=cox_result$coefficients[,2], lower_ci=cox_result$conf.int[,3], upper_ci=cox_result$conf.int[,4], P=cox_result$coefficients[,5])




### ED 5D
Relapse_survdata <- left_join(Relapse_survdata, all_tumour_df_20220209[,c("tumour_id_mphase_cruk","age","sex","pack_years_calculated","pTNMStage_v8","adjuvant_treatment_YN")], by = "tumour_id_mphase_cruk")

#factor levels
Relapse_survdata$sex <- factor(Relapse_survdata$sex , levels = c("Male","Female"))
Relapse_survdata$adjuvant_treatment_YN <- factor(Relapse_survdata$adjuvant_treatment_YN , levels = c("No adjuvant","Adjuvant"))
Relapse_survdata$TNM_combine <- NA
Relapse_survdata$TNM_combine[which(grepl("1",Relapse_survdata$pTNMStage_v8 ))] <- "I"
Relapse_survdata$TNM_combine[which(grepl("2",Relapse_survdata$pTNMStage_v8 ))] <- "II"
Relapse_survdata$TNM_combine[which(grepl("3",Relapse_survdata$pTNMStage_v8 ))] <- "III"
Relapse_survdata$TNM_combine <- factor(Relapse_survdata$TNM_combine, levels = c("I","II","III"))

#MVA
MVA.cox.res <- coxph(formula = Surv(Relapse_survdata$dfs_time, Relapse_survdata$cens_dfs) ~ 
                       Relapse_survdata$ORACLE_mean + 
                       Relapse_survdata$sex + 
                       Relapse_survdata$age + 
                       Relapse_survdata$pack_years_calculated + 
                       Relapse_survdata$adjuvant_treatment_YN +
                       Relapse_survdata$TNM_combine )

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

#
MVA.cox.res$pvalue <-signif(MVA.cox.res$P,digits = 1) 
MVA.cox.res$pvalue[which(is.na(MVA.cox.res$pvalue))] <- "Ref"
MVA.cox.res$mean <- round(MVA.cox.res$HR,digits = 2)

for(i in 1:nrow(MVA.cox.res)){
  lci <- round(MVA.cox.res$lower_ci[i],digits=2);
  uci <- round(MVA.cox.res$upper_ci[i],digits=2);
  MVA.cox.res$result[i] <- paste(lci,uci,sep = " - ")
}

#
MVA.cox.res$mean[which(MVA.cox.res$pvalue == "Ref")] <- ""

#
MVA.cox.res <- rbind(c(1,NA,NA,"","","p-value","HR","CI"),MVA.cox.res)

#turn data from chr to numeric
MVA.cox.res[,1:4] <- as.data.frame(apply(MVA.cox.res[,1:4],2,as.numeric))

#Forest plot
forestplot(MVA.cox.res[,c(5:8)],graph.pos=2,
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





### 4A

#bin the correlation with ORACLE
IC50_results$cor <- ifelse(IC50_results$estimate > 0 , "Positive","Negative")

#label FDA approved drugs for NSCLC
IC50_results$label <- IC50_results$cor
IC50_results$label[which(IC50_results$significance=="FALSE")] <- NA
IC50_results$label[which(IC50_results$FDA=="TRUE")] <- "FDA"

#volcano plot of drug sensitivity to ORACLE risk-score
ggplot(IC50_results,aes(x = estimate,y = -log10(pvalue)))+
  geom_point(alpha=0.7,aes(col=significance),size=3)+
  geom_point(data = IC50_results[which(IC50_results$FDA=="TRUE"),],alpha=0.8,size=3,aes(x = estimate,y = -log10(pvalue)),pch=21,col="black",fill="#D1D1D1")+
  geom_point(data = IC50_results[which(IC50_results$Name=="cisplatin"),],alpha=0.8,size=3,aes(x = estimate,y = -log10(pvalue)),pch=21,col="black",fill="#E20000")+
  scale_x_continuous(breaks=seq(-0.5,0.5,0.5),expand=c(0,0),limits = c(-0.5,0.5))+
  ylab("-log(p-value)") + xlab("") +
  ggtitle("")+
  scale_y_continuous(expand = c(0,0),limits = c(0,3))+
  geom_text_repel(data = IC50_results[which(IC50_results$FDA=="TRUE"),],aes(label = Name),min.segment.length = unit(0.1, "lines"), size=3,max.overlaps = 20)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = -log10(0.05),lty="dashed") + 
  scale_color_manual(values=c("#D1D1D1","#E20000"),na.value = "#D1D1D1")+
  theme(legend.text = element_text(size=12),axis.title = element_text(size=14),axis.text = element_text(size = 12),panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())





### 4B
#correlation between IC50 and ORACLE risk-score classified by targeting pathway
ggplot(IC50_results, aes(fct_reorder(pathway,estimate),estimate)) +  
  geom_point(aes(col=significance),size=3,alpha=0.8) + 
  geom_boxplot(width=0.4,col="#0060A2",alpha=0,outlier.colour = NULL) +
  scale_color_manual(values = c("azure4","#E20000")) + 
  scale_y_continuous(expand = c(0,0),limits = c(-0.5,0.5)) +
  ylab("") + xlab("Targeting pathway") + 
  geom_hline(yintercept = 0, lty="dashed") +
  theme_bw() +
  theme(legend.position = "none",axis.text.x = element_text(angle=90, vjust=0.5,hjust = 1),panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())





### ED 6B

#GDSC IC50 data obtained from DepMap
LUAD_IC50 <- read.csv("LUAD_IC50_20220112.csv",row.names = 1)
df <- read_excel("LUAD_IC50_20220112.xlsx")
colnames(LUAD_IC50) <- colnames(df)[-1]

#plot correlation between IC50 and ORACLE risk-score for significantly-associated drugs
tmp <- IC50_results[which(IC50_results$significance == "TRUE"),]
sig.drugs <- fct_reorder(rownames(tmp), -tmp$estimate)
sig.drugs <- levels(sig.drugs)

for(i in 1:length(sig.drugs)){
  md <- cor.test(LUAD_IC50[,sig.drugs[i]],LUAD_IC50$riskscore,method="spearman",exact=FALSE)
  volc_p <- ggplot(LUAD_IC50, aes(LUAD_IC50[,sig.drugs[i]], riskscore)) + geom_point(size=5) + geom_smooth(method="lm",se=F) +
    xlab("Log(IC50)") + ylab("Risk score") + ggtitle(sig.drugs[i]) +
    #scale_x_continuous(breaks = seq(0,8,2),expand = c(0,0),limits = c(0,9)) +
    scale_y_continuous(expand = c(0,0),limits = c(10.5,13.5)) +
    geom_text(aes(x=min(LUAD_IC50[,sig.drugs[i]], na.rm=T) + 1.4, y=13.4, label = paste0("R=",signif(md$estimate,digits = 2),", P=",signif(md$p.value,digits = 2)))) +
    theme(panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())
  
  print(volc_p)
}

#output significant drugs' info
sig.drugs <- data.frame(drugs = sig.drugs)
write.csv(sig.drugs, "ORACLE_sensitive_drugs.csv", row.names = F)







##load data
Tx421_riskscore <- read.csv("20230727_TRACERx421_riskscores.csv")

#cut-off developed in 2019 Biswas et al.
cutoff <- read.csv("2018-11-07_oracle_de.novo_cut.off.csv")
cutoff <- cutoff$RiskScore[1]

#separate region info
Tx421_riskscore$RegionalID <- apply(Tx421_riskscore, 1, FUN = function(x){ strsplit(x[1], split = "-")[[1]][2] })
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

#
TRACERx_clinical <- readRDS("all_patient_df_20221109.RDS")

all_tumour_df_20220209 <- readRDS('all_tumour_df_20220209.RDS')
all_tumour_df_20220209$tumour_id_mphase_cruk <- apply(all_tumour_df_20220209, 1, FUN = function(x){ if(grepl("Cluster",x[12])) { gsub(".*-", x[12], replacement= paste0(x[4],"_")) }
  else{ x[12] <- x[4] } })





### 5A

exploratory_RiskScore <- Tx421_riskscore

#create df with ORACLE mean RS
exploratory_RiskScore$mean <- NA
for(i in 1:nrow(exploratory_RiskScore)){
  exploratory_RiskScore$mean[i] <- exploratory_RiskScore$ORACLE_riskscore[which(exploratory_RiskScore$cruk_id == exploratory_RiskScore$cruk_id[i])] %>% mean()
}

#retrieve clinical factors
Features.df <- left_join(exploratory_RiskScore,all_tumour_df_20220209[,c("tumour_id_mphase_cruk","patient_3pad","age","sex","smoking_status","pTNMStage_v8")],by="tumour_id_mphase_cruk")
Features.df$pTNMStage_v8[which(Features.df$cruk_id == "CRUK0704")] <- "2b"
Features.df <- Features.df[-which(duplicated(Features.df$cruk_id)),] #remove duplicates

#Biopsy number
Features.df$biopsy_no <- NA
for(i in 1:nrow(Features.df)){
  Features.df$biopsy_no[i] <- sum(Features.df$High[i],Features.df$Low[i])
}

#Smoke status
Features.df$smoking_status[which(grepl("Recent",Features.df$smoking_status))] <- "Ex-Smoker"
Features.df$smoking_status <- factor(Features.df$smoking_status,levels = c("Never Smoked","Ex-Smoker","Current Smoker"))

#Volume
TRACERx_Volume_Pyradiomics <- read.csv("TRACERx_Volume_Pyradiomics.csv")
TRACERx_Volume_Pyradiomics$patient_3pad <- gsub("LTX0",TRACERx_Volume_Pyradiomics$patient_3pad,replacement = "LTX")
TRACERx_Volume_Pyradiomics <- left_join(TRACERx_Volume_Pyradiomics, all_tumour_df_20220209[which(!duplicated(all_tumour_df_20220209$patient_3pad)), c("cruk_id","patient_3pad")], by = "patient_3pad")
TRACERx_Volume_Pyradiomics$tumour_id_mphase_cruk <- apply(TRACERx_Volume_Pyradiomics, 1, FUN = function(x){ if(grepl("Cluster",x[2])){ gsub(".*_", x[2], replacement = paste0( x[24],"_")) } else { x[2] <- x[24] }})
TRACERx_Volume_Pyradiomics$Volume <- TRACERx_Volume_Pyradiomics$Volume/1000
Features.df <- left_join(Features.df,TRACERx_Volume_Pyradiomics[,c("tumour_id_mphase_cruk","Volume")], by = "tumour_id_mphase_cruk")

#KI-67
KI67 <- filter(dplyr::select(all_tumour_df_20220209, cruk_id, KI67_score), !duplicated(cruk_id))
Features.df <- left_join(Features.df,KI67[,c("cruk_id","KI67_score")], by = "cruk_id")

#combine Stage
Features.df$pTNMStage_v8[which(grepl("1",Features.df$pTNMStage_v8))] <- "I"
Features.df$pTNMStage_v8[which(grepl("2",Features.df$pTNMStage_v8))] <- "II"
Features.df$pTNMStage_v8[which(grepl("3",Features.df$pTNMStage_v8))] <- "III"

#prepare evo metrics
Tx421_evo_metrics <- read.table("20221110_TRACERx421_evolutionary_metrics.tsv", header = T, sep = '\t')
Tx421_evo_metrics$tumour_id <- gsub("_Tumour", Tx421_evo_metrics$tumour_id, replacement = "_Cluster")
Tx421_evo_metrics <- Tx421_evo_metrics[which(Tx421_evo_metrics$tumour_id %in% exploratory_RiskScore$tumour_id_mphase_cruk),]
Tx421_evo_metrics$cruk_id <- gsub("_.*",Tx421_evo_metrics$tumour_id,replacement = "")
Tx421_evo_metrics_select <- dplyr::select(Tx421_evo_metrics,"tumour_id","cruk_id",
                                          "wFLOH","frac_abberant_genom_subcl",
                                          "num_clonal_gds","num_subclonal_gds")

#load mutation data
count_drivers <- read.csv("conut_drivers.csv")
count_drivers <- dplyr::select(count_drivers, cruk_id,num_clonal_drivers,num_subclonal_drivers)

#load subclonalExpansionScore
diversity_df <- read.csv("20221102_TRACERx421_tumourDivCorrectedtree.csv")
colnames(diversity_df)[1] <- "cruk_id"
diversity_df$tumour_id <- gsub("_Tumour", diversity_df$tumour_id, replacement = "_Cluster")
diversity_df <- dplyr::select(diversity_df,cruk_id, tumour_id,maxCCF_terminalonly_tum_meanccf  )
diversity_df <- diversity_df  %>%
  group_by(tumour_id) %>%  
  mutate(RecentSubclExpansionScore = max(maxCCF_terminalonly_tum_meanccf, na.rm = T)) %>%
  ungroup() %>%
  dplyr::select(cruk_id, RecentSubclExpansionScore)
diversity_df <- diversity_df[which(!duplicated(diversity_df$cruk_id)),]

#combine mutation data and recent subclonal expansion score
tmp <- left_join(Tx421_evo_metrics_select,count_drivers, by = "cruk_id")
tmp <- left_join(tmp, diversity_df, by="cruk_id")

#average evo metric for patients with multiple primary tumours
tmp <- aggregate(tmp[,-c(1:2)], by = list(cruk_id = tmp$cruk_id), FUN = function(x) { median(x, na.rm = T) })

#percentages of FLOH and SCNA-ITH
tmp$wFLOH <- tmp$wFLOH*100
tmp$frac_abberant_genom_subcl <- tmp$frac_abberant_genom_subcl*100

#level for discrete variables
Features.df <- left_join(Features.df,tmp, by ="cruk_id")
Features.df$sex <- factor(Features.df$sex, levels = c("Male","Female"))
Features.df$pTNMStage_v8 <- factor(Features.df$pTNMStage_v8,levels = c("I","II","III"))

#standardisation
Clinical.scale <- dplyr::select(Features.df, cruk_id, mean, age, biopsy_no, Volume, KI67_score)
dummy <- as.data.frame(model.matrix(data = Features.df, ~sex + pTNMStage_v8 + smoking_status))[,-1]
dummy$cruk_id <- Clinical.scale$cruk_id
Clinical.scale <- left_join(Clinical.scale, dummy, by = "cruk_id")
Clinical.scale[,-1] <- scale(Clinical.scale[,-1], center = T)

Genetics.scale <-  dplyr::select(Features.df, cruk_id, mean, frac_abberant_genom_subcl, wFLOH, num_clonal_gds, num_subclonal_gds, num_clonal_drivers, num_subclonal_drivers, RecentSubclExpansionScore)
Genetics.scale[,-1] <- scale(Genetics.scale[,-1], center = T)


#Multiple linear regression - clinical features
MLR_model_clinical <- lm(data=Clinical.scale,mean~age+sexFemale+`smoking_statusEx-Smoker`+`smoking_statusCurrent Smoker`+ `pTNMStage_v8II`+ `pTNMStage_v8III` + biopsy_no+Volume+KI67_score)
r <- summary(MLR_model_clinical)
pvalue.MLR.clinical <- data.frame(correlate = names(r$coefficients[-1,1]),pvalue = r$coefficients[-1,4],coef = r$coefficients[-1,1])
pvalue.MLR.clinical$correlate <- c("Age","Female","Ex-Smoker","Smoker","Stage II","Stage III","#biopsy","Volume","Ki67")
pvalue.MLR.clinical$correlate <- factor(pvalue.MLR.clinical$correlate , levels = pvalue.MLR.clinical$correlate)

pvalue.MLR.clinical$significance <- ""
pvalue.MLR.clinical$significance[which(pvalue.MLR.clinical$pvalue < 0.05 & pvalue.MLR.clinical$pvalue >= 0.01)] <- "*"
pvalue.MLR.clinical$significance[which(pvalue.MLR.clinical$pvalue < 0.01 & pvalue.MLR.clinical$pvalue >= 0.005)] <- "**"
pvalue.MLR.clinical$significance[which(pvalue.MLR.clinical$pvalue < 0.005)] <- "***"

mean_clin <- ggplot(pvalue.MLR.clinical)+
  geom_tile(aes(correlate, y = 1,fill=coef))+
  scale_y_continuous(expand = c(0,0))+
  xlab("")+ylab("Mean")+labs(col="")+
  scale_fill_distiller(name = '',limits=c(-0.45,0.45),palette ="RdBu",na.value = "white") +
  scale_x_discrete(expand = c(0,0)) +
  geom_text(aes(x = correlate, y = 1, label=significance),size=5,hjust=0.5)+
  theme(legend.position = "top",axis.text = element_blank(), axis.title.y = element_text(angle = 0, vjust=0.5,hjust = 1),axis.ticks = element_blank(),panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())


#Multiple linear regression - genetics features
MLR_model_genetic <- lm(data = Genetics.scale, mean ~  frac_abberant_genom_subcl + wFLOH + num_clonal_gds + num_subclonal_gds + num_clonal_drivers + num_subclonal_drivers + RecentSubclExpansionScore  ) 

r <- summary(MLR_model_genetic)
pvalue.MLR.genetics <- data.frame(correlate = names(r$coefficients[-1,1]),pvalue = r$coefficients[-1,4],coef = r$coefficients[-1,1])
pvalue.MLR.genetics$correlate <- c("SCNA-ITH","FLOH","Clonal WGD","Subclonal WGD","Clonal drivers","Subclonal drivers","Recent subclonal expansion")
pvalue.MLR.genetics$correlate <- factor(pvalue.MLR.genetics$correlate , levels = pvalue.MLR.genetics$correlate)

pvalue.MLR.genetics$significance <- ""
pvalue.MLR.genetics$significance[which(pvalue.MLR.genetics$pvalue < 0.1 & pvalue.MLR.genetics$pvalue >= 0.01)] <- "*"
pvalue.MLR.genetics$significance[which(pvalue.MLR.genetics$pvalue < 0.01 & pvalue.MLR.genetics$pvalue >= 0.005)] <- "**"
pvalue.MLR.genetics$significance[which(pvalue.MLR.genetics$pvalue < 0.005)] <- "***"

mean_gene <- ggplot(pvalue.MLR.genetics)+
  geom_tile(aes(correlate, y = 1,fill=coef))+
  scale_y_continuous(expand = c(0,0))+
  xlab("")+ylab("")+labs(col="")+
  scale_fill_distiller(name = '',limits=c(-0.45,0.45),palette ="RdBu",na.value = "white") +
  scale_x_discrete(expand = c(0,0)) +
  geom_text(aes(x = correlate, y = 1, label=significance),size=5,hjust=0.5)+
  theme(legend.position = "top",axis.text = element_blank(), axis.title.y = element_text(angle = 0, vjust=0.5,hjust = 1),axis.ticks = element_blank(),panel.border = element_rect(colour = "black",fill = NA,size=1),panel.background = element_blank())



### ED 8A
#plot correlation between ORACLE risk-score and clinicopathological/genetic features
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






### 5B

#load all TRACERx biomarker data
Tx_biomarkers <- readRDS("20221109_LUAD_input_plot_df_all.248.RDS")
ctDNA.df <- read.csv("20221108_Tx421_ctDNA_data_LUAD_single.csv")
Evo_metrics <- read.table("20221110_TRACERx421_evolutionary_metrics.tsv", header = T, sep = '\t')

#Identify STAS status per patient tumour
Tx_biomarkers <- dplyr::select(Tx_biomarkers, tumour_id_mphase_cruk,cruk_id,STAS_presence) 
Tx_biomarkers$tumour_id_mphase_cruk <- gsub("_Tumour",Tx_biomarkers$tumour_id_mphase_cruk,replacement = "_Cluster")
Tx_biomarkers <- Tx_biomarkers[which(Tx_biomarkers$tumour_id_mphase_cruk %in% Tx421_riskscore$tumour_id_mphase_cruk),]
Tx_biomarkers <- Tx_biomarkers %>% group_by(cruk_id) %>%
  summarise(STAS = case_when(any(STAS_presence == "present") ~ "present",
                             all(STAS_presence == "absent") ~ "absent")
  ) %>% ungroup()

Tx_biomarkers <- Tx_biomarkers[which(!duplicated(Tx_biomarkers$cruk_id)),]

#ctDNA status
ctDNA.df <- ctDNA.df %>% mutate(ctDNA_status = case_when(preop_ctDNA_patient == "yes" | preop_ctDNA_natera == "pos" ~ "shedder",
                                                         preop_ctDNA_patient == "no" | preop_ctDNA_natera == "neg" ~ "low-shedder"))

#obtain evo metrics
Evo_metrics$tumour_id <- gsub("_Tumour", Evo_metrics$tumour_id, replacement = "_Cluster")
Evo_metrics <- Evo_metrics[which(Evo_metrics$tumour_id %in% Tx421_riskscore$tumour_id_mphase_cruk),]
Evo_metrics$cruk_id <- gsub("_.*",Evo_metrics$tumour_id,replacement = "")
Evo_metrics <- dplyr::select(Evo_metrics,"tumour_id","cruk_id",
                             "frac_abberant_genom_subcl","num_clonal_gds","num_subclonal_gds")
Evo_metrics <- Evo_metrics %>% group_by(cruk_id) %>%
  summarise(subclonal_wgd = case_when(num_subclonal_gds == 0 ~ "No",
                                      num_subclonal_gds > 0 ~ "Yes"),
            SCNA_ITH = mean(frac_abberant_genom_subcl, na.rm=T))
Evo_metrics <- Evo_metrics[which(!duplicated(Evo_metrics$cruk_id)),]

#load subclonalExpansionScore
diversity_df <- read.csv("20221102_TRACERx421_tumourDivCorrectedtree.csv")
colnames(diversity_df)[1] <- "cruk_id"
diversity_df$tumour_id <- gsub("_Tumour", diversity_df$tumour_id, replacement = "_Cluster")
diversity_df <- dplyr::select(diversity_df,cruk_id, tumour_id,maxCCF_terminalonly_tum_meanccf  )
diversity_df <- diversity_df  %>%
  group_by(tumour_id) %>%  
  mutate(RecentSubclExpansionScore = max(maxCCF_terminalonly_tum_meanccf, na.rm = T)) %>%
  ungroup() %>%
  dplyr::select(cruk_id, RecentSubclExpansionScore)
diversity_df <- diversity_df[which(!duplicated(diversity_df$cruk_id)),]

#join all biomarkers into one df
Tx_biomarkers <- left_join(Tx_biomarkers, ctDNA.df[,c("cruk_id","ctDNA_status")], by = "cruk_id")
Tx_biomarkers <- left_join(Tx_biomarkers, Evo_metrics,by="cruk_id")
Tx_biomarkers <- left_join(Tx_biomarkers, diversity_df, by = "cruk_id")

#add ORACLE mean risk-score
tmp <- Tx421_riskscore
for(i in 1:nrow(tmp)){
  tmp$ORACLE_mean[i] <- tmp$ORACLE_riskscore[which(tmp$cruk_id == tmp$cruk_id[i])] %>% mean
}
tmp <- tmp[which(!duplicated(tmp$cruk_id)),]
Tx_biomarkers <- left_join(Tx_biomarkers, tmp, by ="cruk_id")

#obtain survival data
Tx_biomarkers <- left_join(Tx_biomarkers, TRACERx_clinical[,c("cruk_id","cens_os","os_time","age","sex","pTNMStage_v8","adjuvant_treatment_YN","pack_years_calculated")], by ="cruk_id")
Tx_biomarkers$os_time <- Tx_biomarkers$os_time/365
Tx_biomarkers$cens_os[which(Tx_biomarkers$os_time >5)]<-0
Tx_biomarkers$os_time[which(Tx_biomarkers$os_time>5)] <- 5

#combine stages
Tx_biomarkers$pTNMStage_v8<- as.character(Tx_biomarkers$pTNMStage_v8)
Tx_biomarkers$pTNMStage_v8[which(grepl("1",Tx_biomarkers$pTNMStage_v8))] <- "I"
Tx_biomarkers$pTNMStage_v8[which(grepl("2",Tx_biomarkers$pTNMStage_v8))] <- "II"
Tx_biomarkers$pTNMStage_v8[which(grepl("3",Tx_biomarkers$pTNMStage_v8))] <- "III"

#remove patients with any data point unavailable
Tx_biomarker_avail <- Tx_biomarkers[which(!is.na(Tx_biomarkers$STAS)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$ctDNA_status)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$SCNA_ITH)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$subclonal_wgd)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$RecentSubclExpansionScore)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$ORACLE_mean)),]

#create dummy variables for biomarker measurements
Tx_biomarker_avail$STAS[which(Tx_biomarker_avail$STAS == "present")] <- 1
Tx_biomarker_avail$STAS[which(Tx_biomarker_avail$STAS == "absent")] <- 0
Tx_biomarker_avail$ctDNA_status[which(Tx_biomarker_avail$ctDNA_status == "low-shedder")] <- 0
Tx_biomarker_avail$ctDNA_status[which(Tx_biomarker_avail$ctDNA_status == "shedder")] <- 1
Tx_biomarker_avail$subclonal_wgd[which(Tx_biomarker_avail$subclonal_wgd == "No")] <- 0
Tx_biomarker_avail$subclonal_wgd[which(Tx_biomarker_avail$subclonal_wgd == "Yes")] <- 1
Tx_biomarker_avail$adjuvant_treatment_YN[which(Tx_biomarker_avail$adjuvant_treatment_YN == "No adjuvant")] <- 0
Tx_biomarker_avail$adjuvant_treatment_YN[which(Tx_biomarker_avail$adjuvant_treatment_YN == "Adjuvant")] <- 1

#level
Tx_biomarker_avail$STAS <- factor(Tx_biomarker_avail$STAS, levels = c(0,1))
Tx_biomarker_avail$subclonal_wgd <- factor(Tx_biomarker_avail$subclonal_wgd, levels = c(0,1))
Tx_biomarker_avail$ctDNA_status <- factor(Tx_biomarker_avail$ctDNA_status, levels = c(0,1))
Tx_biomarker_avail$adjuvant_treatment_YN <- factor(Tx_biomarker_avail$adjuvant_treatment_YN, levels = c("No adjuvant","Adjuvant"))
Tx_biomarker_avail$pTNMStage_v8 <- factor(Tx_biomarker_avail$pTNMStage_v8, levels = c("I","II","III"))
Tx_biomarker_avail$sex <- factor(Tx_biomarker_avail$sex, levels = c("Male","Female"))


#MVA
MVA.cox.res <- coxph(data = Tx_biomarker_avail,formula = Surv(Tx_biomarker_avail$os_time,Tx_biomarker_avail$cens_os) ~ 
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

ggforest(data = as.data.frame(Tx_biomarker_avail),MVA.cox.res)





### 5C

#remove patients with any data point unavailable
Tx_biomarker_avail <- Tx_biomarkers[which(!is.na(Tx_biomarkers$STAS)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$ctDNA_status)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$SCNA_ITH)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$subclonal_wgd)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$RecentSubclExpansionScore)),]
Tx_biomarker_avail <- Tx_biomarker_avail[which(!is.na(Tx_biomarker_avail$ORACLE_mean)),]


#binary code
Tx_biomarker_avail$STAS[which(Tx_biomarker_avail$STAS == "absent")] <- 0
Tx_biomarker_avail$STAS[which(Tx_biomarker_avail$STAS == "present")] <- 1

Tx_biomarker_avail$ctDNA_status[which(Tx_biomarker_avail$ctDNA_status == "low-shedder")] <- 0
Tx_biomarker_avail$ctDNA_status[which(Tx_biomarker_avail$ctDNA_status == "shedder")] <- 1

Tx_biomarker_avail$subclonal_wgd[which(Tx_biomarker_avail$subclonal_wgd == "No")] <- 0
Tx_biomarker_avail$subclonal_wgd[which(Tx_biomarker_avail$subclonal_wgd == "Yes")] <- 1

#scale ORACLE score
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
Tx_biomarker_avail$ORACLE_scaled <- range01(Tx_biomarker_avail$ORACLE_mean)

#perform glm test
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

##iterate from FU time = 1 to 5
glm_time_all <- list()
for(i in 1:5){
  glm_time_all[[i]] <- biomarker.glm(i)
}

glm_time_all <- Reduce(rbind,glm_time_all)

#bar plot
ggplot(glm_time_all, aes(OS_time, Explained_var)) + xlab("Time (years)") + ylab("% variance explained") + 
  geom_col(aes(fill=variable),col="black",alpha=0.8) +
  scale_fill_manual(values = brewer.pal(12,"Paired")[c(7,1,4,9,6,3)]) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,40)) +
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))






### ED 8B

#load mutation data
Tx421_regionMut <- fst::read.fst("tx421_RegionMutTable.fst")

#filter for driver mutations
Tx421_regionMut <- Tx421_regionMut[which(Tx421_regionMut$PASS == "TRUE" & Tx421_regionMut$Is.present.region == "TRUE"),]
Tx421_regionMut <- Tx421_regionMut[which(Tx421_regionMut$DriverMut == "TRUE"),]

#fix names
Tx421_regionMut$RegionID <- gsub("_LN0",Tx421_regionMut$RegionID,replacement = "_LN")
Tx421_regionMut$RegionID <- gsub("_LN",Tx421_regionMut$RegionID,replacement = "_LN0")
Tx421_regionMut$RegionID <- gsub(":SU_T1.",Tx421_regionMut$RegionID,replacement = "_SU_T1-")
Tx421_regionMut$RegionID <- gsub(":SU_T2.",Tx421_regionMut$RegionID,replacement = "_SU_T2-")
Tx421_regionMut$RegionID <- gsub(":SU_LN",Tx421_regionMut$RegionID,replacement = "_SU_LN")

#clonal
Clonal.SNV <- Tx421_regionMut[which(Tx421_regionMut$PyCloneClonal_SC == "C"),c("RegionID","PyCloneClonal_SC","Hugo_Symbol")]
Clonal.SNV <- as.data.frame(dcast(Clonal.SNV, RegionID~Hugo_Symbol, value.var = "PyCloneClonal_SC"))
colnames(Clonal.SNV)[1] <- "sample_name"

#binary df for all primary Tx samples
clonal.freq.df <- data.frame(sample_name = gsub("LTX0",Tx421_riskscore$sample_name,replacement = "LTX"), bin = Tx421_riskscore$bin)
clonal.freq.df <- left_join(clonal.freq.df, Clonal.SNV)
for(i in 3:ncol(clonal.freq.df)){
  clonal.freq.df[which(is.na(clonal.freq.df[,i])),i] <- 0
}

for(i in 3:ncol(clonal.freq.df)){
  clonal.freq.df[which(clonal.freq.df[,i] > 1),i] <- 1
}

#
rownames(clonal.freq.df) <- clonal.freq.df$sample_name
clonal.freq.df <- dplyr::select(clonal.freq.df, -contains("sample_name"))

#remove low frequency genes (<2%)
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

#annotate the significance
clonal.result.df$significance <- ifelse(clonal.result.df$pvalue<0.05,"Yes","No")
clonal.result.df$enrichment <- ifelse(clonal.result.df$odds_ratio < 1,"Concordant low\nenriched","Concordant high\nenriched")
clonal.result.df$enrichment[which(clonal.result.df$significance=="No")] <- "non-significant"

ggplot(clonal.result.df,aes(log(odds_ratio),-log10(pvalue)))+geom_point(aes(fill=enrichment),pch=21,size=3)+
  scale_x_continuous(expand = c(0,0),limits = c(-2,2),oob = squish)+
  scale_y_continuous(breaks = seq(0,6,1),expand = c(0,0),limits = c(0,6))+
  scale_fill_manual(values = c("#EE0000FF","#3B4992FF","azure4"))+
  ggtitle("Clonal drivers in low vs high")+
  geom_text_repel(data = subset(clonal.result.df,significance=="Yes"),aes(label=Driver),max.overlaps = 25,size=3)+
  labs(fill="")+xlab("Odds ratio")+ylab("-log10(p-value)")+
  geom_vline(xintercept = 0)+geom_hline(yintercept = -log10(0.05),lty="dashed")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))



#subclonal mutation data
Subclonal.SNV <- Tx421_regionMut[which(Tx421_regionMut$PyCloneClonal_SC == "S"),c("RegionID","PyCloneClonal_SC","Hugo_Symbol")]
Subclonal.SNV <- as.data.frame(dcast(Subclonal.SNV, RegionID~Hugo_Symbol, value.var = "PyCloneClonal_SC"))
colnames(Subclonal.SNV)[1] <- "sample_name"

#binary df for all primary Tx samples
Subclonal.freq.df <- data.frame(sample_name = gsub("LTX0",Tx421_riskscore$sample_name,replacement = "LTX"), bin = Tx421_riskscore$bin)
Subclonal.freq.df <- left_join(Subclonal.freq.df, Subclonal.SNV)
for(i in 3:ncol(Subclonal.freq.df)){
  Subclonal.freq.df[which(is.na(Subclonal.freq.df[,i])),i] <- 0
}

for(i in 3:ncol(Subclonal.freq.df)){
  Subclonal.freq.df[which(Subclonal.freq.df[,i] > 1),i] <- 1
}

#
rownames(Subclonal.freq.df) <- Subclonal.freq.df$sample_name
Subclonal.freq.df <- dplyr::select(Subclonal.freq.df, -contains("sample_name"))

#remove low frequency genes (<1%)
Subclonal.freq.df <- Subclonal.freq.df[, -c(which(colSums(Subclonal.freq.df[,-1]) == 0) + 1 )]
Subclonal.freq.df <- Subclonal.freq.df[,-c(which(colSums(Subclonal.freq.df[,-1]) < nrow(Subclonal.freq.df)*0.01)+1)]
Subclonal.freq.df$bin <- factor(Subclonal.freq.df$bin, levels = c("Low","High"))

#Fisher exact test pvalue
Subclonal.result.df <- data.frame(Driver = colnames(Subclonal.freq.df)[-1],pvalue=NA,odds_ratio = NA)

for(i in 1:nrow(Subclonal.result.df)){
  Subclonal.result.df$pvalue[i] <- fisher.test(table(Subclonal.freq.df$bin,Subclonal.freq.df[,i+1]))$p.value
  Subclonal.result.df$odds_ratio[i] <- fisher.test(table(Subclonal.freq.df$bin,Subclonal.freq.df[,i+1]))$estimate
}


#annotate the significance
Subclonal.result.df$significance <- ifelse(Subclonal.result.df$pvalue<0.05,"Yes","No")
Subclonal.result.df$enrichment <- ifelse(Subclonal.result.df$odds_ratio < 1,"Concordant low\nenriched","Concordant high\nenriched")
Subclonal.result.df$enrichment[which(Subclonal.result.df$significance=="No")] <- "non-significant"

ggplot(Subclonal.result.df,aes(log(odds_ratio),-log10(pvalue)))+geom_point(aes(fill=enrichment),pch=21,size=3)+
  scale_x_continuous(expand = c(0,0),limits = c(-2,2),oob = squish)+
  scale_y_continuous(breaks = seq(0,6,1),expand = c(0,0),limits = c(0,6))+
  scale_fill_manual(values = c("#EE0000FF","azure4","#3B4992FF"))+
  ggtitle("Subclonal drivers in low vs high")+
  geom_text_repel(data = subset(Subclonal.result.df,significance=="Yes"),aes(label=Driver),max.overlaps = 25,size=3)+
  labs(fill="")+xlab("Odds ratio")+ylab("-log10(p-value)")+
  geom_vline(xintercept = 0)+geom_hline(yintercept = -log10(0.05),lty="dashed")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))






### ED 8C

#load GISTIC results

HighLow_summary <- read.csv("20230620_concordant_oracle_highlow_cytoband_summary_wide.csv")
HighLow_summary$high_sig <- NA
HighLow_summary$high_sig[which(HighLow_summary$log_qval_high > -log10(0.05) & HighLow_summary$Gscore_diff_raw > 0)] <- HighLow_summary$Type[which(HighLow_summary$log_qval_high > -log10(0.05) & HighLow_summary$Gscore_diff_raw > 0)]
HighLow_summary$high_sig[which(is.na(HighLow_summary$high_sig))] <- "Not significant"

#
HighLow_summary$low_sig <- NA
HighLow_summary$low_sig[which(HighLow_summary$log_qval_low > -log10(0.05) & HighLow_summary$Gscore_diff_raw < 0)] <- HighLow_summary$Type[which(HighLow_summary$log_qval_low > -log10(0.05) & HighLow_summary$Gscore_diff_raw < 0)]
HighLow_summary$low_sig[which(is.na(HighLow_summary$low_sig))] <- "Not significant"

#significantly detected cytobands and chr amrs
high.sig.arm <- HighLow_summary$cytoband_arm[which( HighLow_summary$high_sig != "Not significant")]%>%unique()
high_sig <- HighLow_summary$cytoband_name[which(HighLow_summary$high_sig %in% c("Amp","Del"))] %>% unique()

low.sig.arm <- HighLow_summary$cytoband_arm[which( HighLow_summary$low_sig != "Not significant")]%>%unique()
low_sig <- HighLow_summary$cytoband_name[which(HighLow_summary$low_sig %in% c("Amp","Del"))] %>% unique()

#
annot_genes <- fst::read_fst( "gene_anno_df.fst", as.data.table=T )
cn_drivers <- annot_genes[ (is_cn_driver), gene_name ]
annot_genes$chr_arm <- gsub("p.*", annot_genes$cytoband_name, replacement = "p")
annot_genes$chr_arm <- gsub("q.*", annot_genes$chr_arm, replacement = "q")


#
mphase <- read.fst("20221109_TRACERx421_scna_table.fst")

#segment name
for(i in 1:nrow(mphase)){
  mphase$seg_name[i] <- paste(mphase$chr[i], mphase$startpos[i], mphase$endpos[i], sep = ":")
}
mphase <- mphase[,c("sample","chr","startpos","endpos","cpn_event_vs_ploidy","seg_name","patient_id","tumour_id")]

#get genome position and correspond gene symbols
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
attributes <- c("hgnc_symbol",'chromosome_name', 
                'start_position', 'end_position')
gene_positions <- as.data.table( getBM(attributes=attributes,mart=ensembl, 
                                       useCache = FALSE) )

# filter for annotated nuclear genes 1-22
gene_positions <- gene_positions[ chromosome_name  %in% 1:22 ]

#annotate onto GISTIC output
chosen_gene_positions <- makeGRangesFromDataFrame( gene_positions[ hgnc_symbol %in% cn_drivers, .(chr = chromosome_name, 
                                                                                                  start = start_position, 
                                                                                                  end = end_position, 
                                                                                                  hgnc_symbol) ], keep.extra.columns = TRUE )
#find overlapping genome position with our mphase data
colnames(mphase)[2:4] <- c("chromosome","start","end")
mphase_gr <- makeGRangesFromDataFrame( mphase, keep.extra.columns = TRUE )
overlaps <- findOverlaps(mphase_gr, chosen_gene_positions)
mphase_matched <- cbind( as.data.table(mphase_gr)[ overlaps@from ], 
                         as.data.table(chosen_gene_positions)[ overlaps@to ] )

mphase <- left_join(mphase, mphase_matched[which(!duplicated(mphase_matched$seg_name)),c("seg_name","hgnc_symbol")], by = "seg_name")
colnames(mphase)[7:8] <- c("cruk_id","tumour_id_mphase_cruk")

#select LUAD patients in Tx421 cohort
mphase <- mphase[which(mphase$tumour_id_mphase_cruk %in% Tx421_riskscore$tumour_id_mphase_cruk),]

#calculate frequencies of clonal or subclonal events with segment loss/gain
mphase_freq <- mphase %>% group_by(seg_name,cruk_id) %>%
  summarise(total = length(sample),
            loss.regions = sum(cpn_event_vs_ploidy %in% c('loss','deep_loss')),
            neutral.regions = sum(cpn_event_vs_ploidy == 'neutral'),
            gain.regions = sum(cpn_event_vs_ploidy %in% c('gain','amp'))) %>%
  mutate(loss.subclonal = ifelse(loss.regions < total & loss.regions != 0, TRUE, FALSE)) %>%
  mutate(gain.subclonal = ifelse(gain.regions < total & gain.regions != 0, TRUE, FALSE)) %>%
  mutate(loss.clonal = ifelse(loss.regions == total & loss.regions!= 0 ,TRUE, FALSE)) %>%
  mutate(gain.clonal = ifelse(gain.regions == total & gain.regions!= 0 ,TRUE, FALSE))

#join driver genes and ORACLE risk class
mphase_freq <- left_join(mphase_freq, mphase[which(!duplicated(mphase$seg_name)),c("seg_name","hgnc_symbol")], by = "seg_name")
colnames(mphase_freq)[11] <- "driver"
mphase_freq <- left_join(mphase_freq, Tx421_riskscore[which(!duplicated(Tx421_riskscore$cruk_id)),c("cruk_id","ORACLE_class")], by = "cruk_id")


#annotate arm level SCNA
for(i in 1:nrow(HighLow_summary)){
  if(HighLow_summary$high_sig[i] %in% c("Amp","Del") | HighLow_summary$low_sig[i] %in% c("Amp","Del")){
    HighLow_summary$SCNA[i] <- paste(HighLow_summary$cytoband_arm[i], HighLow_summary$Type[i], sep = ";")
  }
  else{HighLow_summary$SCNA[i] <- "Not significant"}
}

#create df for significant chromosome arm, SCNA type and driver gene located in
High.sig.df <- data.frame(SCNA = unique(HighLow_summary$SCNA[which(HighLow_summary$high_sig != "Not significant")]), Driver = NA, cohort = "High-risk")
High.sig.df$chr_arm <- sapply(High.sig.df$SCNA,FUN = function(x){strsplit(x, split = ";")[[1]][1]})
High.sig.df$Type <- sapply(High.sig.df$SCNA,FUN = function(x){strsplit(x, split = ";")[[1]][2]})

Low.sig.df <- data.frame(SCNA = unique(HighLow_summary$SCNA[which(HighLow_summary$low_sig != "Not significant")]), Driver = NA, cohort = "Low-risk")
Low.sig.df$chr_arm <- sapply(Low.sig.df$SCNA,FUN = function(x){strsplit(x, split = ";")[[1]][1]})
Low.sig.df$Type <- sapply(Low.sig.df$SCNA,FUN = function(x){strsplit(x, split = ";")[[1]][2]})

#fit in driver genes
for(i in 1:nrow(High.sig.df)){
  if(length(which(annot_genes$chr_arm == High.sig.df$chr_arm[i] & annot_genes$is_cn_driver == "TRUE")) > 0){
    High.sig.df$Driver[i] <- paste(paste(unique(annot_genes$gene_name[which(annot_genes$chr_arm == High.sig.df$chr_arm[i] & annot_genes$is_cn_driver == "TRUE")]), High.sig.df$Type[i], sep = "_"), collapse = ";")
  }
  else{ High.sig.df$Driver[i] <- ""}
}

for(i in 1:nrow(Low.sig.df)){
  if(length(which(annot_genes$chr_arm == Low.sig.df$chr_arm[i] & annot_genes$is_cn_driver == "TRUE")) > 0){
    Low.sig.df$Driver[i] <- paste(paste(unique(annot_genes$gene_name[which(annot_genes$chr_arm == Low.sig.df$chr_arm[i] & annot_genes$is_cn_driver == "TRUE")]), Low.sig.df$Type[i], sep = "_"), collapse = ";")
  }
  else{ Low.sig.df$Driver[i] <- ""}
}

#
High.freq.df <- data.frame(driver = rep(unlist(strsplit(High.sig.df$Driver, split = ";")),4),
                           event = rep(c("clonal","subclonal"), each = length(unlist(strsplit(High.sig.df$Driver, split = ";")))), 
                           freq = NA, 
                           cohort = rep(c("High","Low"), each = 16)
)

Low.freq.df <- data.frame(driver = rep(unlist(strsplit(Low.sig.df$Driver, split = ";")),4),
                          event = rep(c("clonal","subclonal"), each = length(unlist(strsplit(Low.sig.df$Driver, split = ";")))), 
                          freq = NA, 
                          cohort = rep(c("High","Low"), each = 50)
)

#
High.freq.df$Type <- sapply(High.freq.df$driver, FUN = function(x){ strsplit(x, split = "_")[[1]][2]})
High.freq.df$driver <- sapply(High.freq.df$driver, FUN = function(x){ strsplit(x, split = "_")[[1]][1]})

Low.freq.df$Type <- sapply(Low.freq.df$driver, FUN = function(x){ strsplit(x, split = "_")[[1]][2]})
Low.freq.df$driver <- sapply(Low.freq.df$driver, FUN = function(x){ strsplit(x, split = "_")[[1]][1]})

#
for(i in 1:nrow(High.freq.df)){
  total_n <- length(which(mphase_freq$driver == High.freq.df$driver[i] & mphase_freq$ORACLE_class == High.freq.df$cohort[i]))
  if(High.freq.df$Type[i] == "Amp" & High.freq.df$event[i] == "clonal"){
    High.freq.df$freq[i] <- length(which(mphase_freq$ORACLE_class == High.freq.df$cohort[i] & mphase_freq$driver == High.freq.df$driver[i] & mphase_freq$gain.clonal == "TRUE"))/total_n
  }
  else if(High.freq.df$Type[i] == "Del" & High.freq.df$event[i] == "clonal"){
    High.freq.df$freq[i] <- length(which(mphase_freq$ORACLE_class == High.freq.df$cohort[i] & mphase_freq$driver == High.freq.df$driver[i] & mphase_freq$loss.clonal == "TRUE"))/total_n
  }
  else if(High.freq.df$Type[i] == "Amp" & High.freq.df$event[i] == "subclonal"){
    High.freq.df$freq[i] <- length(which(mphase_freq$ORACLE_class == High.freq.df$cohort[i] & mphase_freq$driver == High.freq.df$driver[i] & mphase_freq$gain.subclonal == "TRUE"))/total_n
  }
  else if(High.freq.df$Type[i] == "Del" & High.freq.df$event[i] == "subclonal"){
    High.freq.df$freq[i] <- length(which(mphase_freq$ORACLE_class == High.freq.df$cohort[i] & mphase_freq$driver == High.freq.df$driver[i] & mphase_freq$loss.subclonal == "TRUE"))/total_n
  }  
  else{}
}

for(i in 1:nrow(Low.freq.df)){
  total_n <- length(which(mphase_freq$driver == Low.freq.df$driver[i] & mphase_freq$ORACLE_class == Low.freq.df$cohort[i]))
  if(Low.freq.df$Type[i] == "Amp" & Low.freq.df$event[i] == "clonal"){
    Low.freq.df$freq[i] <- length(which(mphase_freq$ORACLE_class == Low.freq.df$cohort[i] & mphase_freq$driver == Low.freq.df$driver[i] & mphase_freq$gain.clonal == "TRUE"))/total_n
  }
  else if(Low.freq.df$Type[i] == "Del" & Low.freq.df$event[i] == "clonal"){
    Low.freq.df$freq[i] <- length(which(mphase_freq$ORACLE_class == Low.freq.df$cohort[i] & mphase_freq$driver == Low.freq.df$driver[i] & mphase_freq$loss.clonal == "TRUE"))/total_n
  }
  else if(Low.freq.df$Type[i] == "Amp" & Low.freq.df$event[i] == "subclonal"){
    Low.freq.df$freq[i] <- length(which(mphase_freq$ORACLE_class == Low.freq.df$cohort[i] & mphase_freq$driver == Low.freq.df$driver[i] & mphase_freq$gain.subclonal == "TRUE"))/total_n
  }
  else if(Low.freq.df$Type[i] == "Del" & Low.freq.df$event[i] == "subclonal"){
    Low.freq.df$freq[i] <- length(which(mphase_freq$ORACLE_class == Low.freq.df$cohort[i] & mphase_freq$driver == Low.freq.df$driver[i] & mphase_freq$loss.subclonal == "TRUE"))/total_n
  }  
  else{}
}

#
High.freq.df$SCNA <- apply(High.freq.df, 1, FUN = function(x){ paste(x[1], x[5], sep = "_")})
Low.freq.df$SCNA <- apply(Low.freq.df, 1, FUN = function(x){ paste(x[1], x[5], sep = "_")})

#
tmp <- Low.freq.df[-which(Low.freq.df$SCNA %in% High.freq.df$SCNA),]
tmp <- rbind(tmp, High.freq.df)
tmp <- tmp %>% group_by(SCNA,cohort) %>% 
  mutate(total_freq = signif(sum(freq), digits = 2))


#genome wide plot - high risk cohort
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
