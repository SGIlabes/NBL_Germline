
# 1. Library ----------------------------------------------------------------------------------
package_list <- c("readr", 'tidyverse',"tibble", "openxlsx", "readxl", 
                  "writexl", "stringr", "naniar", "MASS", "devtools", 
                  "ggvenn", "cowplot", "ggpubr", "survival", "survminer",
                  "hrbrthemes", "viridis", "RColorBrewer", "scales", "psych",'smplot2')
purrr::walk(package_list, library, character.only = TRUE)

output_folder <- '/path/output/TARGET'
data_folder <-'/path/data/TARGET/'

# 2. Data -------------------------------------------------------------------------------------
TARGET_data <- readRDS(paste0(data_folder,'/TARGET_data.RDS'))


# 3. Ethnic -----------------------------------------------------------------------------------
## EthSeq 
# recal.vcf --> hg19tohg38 --> duplication rm (Chromosome, Position, Alternative Allele) d/t liftover --> Merge 


# 4. Germline Functional ----------------------------------------------------------------------

allfilelist <- list.files(path = data_folder, pattern = "target", recursive = FALSE, full.names = TRUE)
filecnt = length(allfilelist)

raw_variant <- data.frame()
for (i in 1:filecnt) {
  # Read the data
  tmpdata <- read.delim(allfilelist[i])
  
  # Rename the 86th column to "GT"
  colnames(tmpdata)[86] <- "GT"
  
  # Split AD column into n_ref_count and t_alt_count
  ad_counts <- str_split_fixed(tmpdata$AD, ',', 2)
  tmpdata$n_ref_count <- as.numeric(ad_counts[,1])
  tmpdata$t_alt_count <- as.numeric(ad_counts[,2])
  
  # Ensure numeric conversion for specific columns
  tmpdata$VAF <- as.numeric(tmpdata$VAF)
  tmpdata$DP <- as.numeric(tmpdata$DP)
  tmpdata$GQ <- as.numeric(tmpdata$GQ)
  tmpdata$AF_gnomad <- as.numeric(tmpdata$AF_gnomad)
  
  # Extract resource ID info from the filename
  resource_id_info <- str_sub(allfilelist[i], 75, 83)
  resource_id_column <- data.frame(germline_id = resource_id_info)
  
  # Apply filters to the data
  tmpdata <- tmpdata %>% 
    filter(DP >= 15) %>% 
    filter(GQ >= 50) %>% 
    filter(VAF >= 0.2) %>% 
    filter(!str_detect(FILTER, "PASS")) %>%  #VQSRTrancheSNP99.9
    filter(FILTER != ".") %>% 
    filter(Chr != "chrY")
  
  # Combine resource ID with the filtered data
  binded1 <- cbind(resource_id_column, tmpdata)
  binded1<-left_join(binded1,TARGET_data[,c("germline_id","Estimated population")],by="germline_id")
  binded1$Estimated_population <-binded1$`Estimated population`
  binded1 <- binded1 %>%
    mutate(ExAC_Filtered = case_when(
      Estimated_population == 'EUR' ~ ifelse(ExAC_ALL < 0.01 | is.na(ExAC_ALL), TRUE, FALSE),
      Estimated_population == 'AFR' ~ ifelse(ExAC_AFR < 0.01 | is.na(ExAC_AFR), TRUE, FALSE),
      Estimated_population == 'EAS' ~ ifelse(ExAC_EAS < 0.01 | is.na(ExAC_EAS), TRUE, FALSE),
      Estimated_population == 'SAS' ~ ifelse(ExAC_SAS < 0.01 | is.na(ExAC_SAS), TRUE, FALSE),
      Estimated_population == 'AMR' ~ ifelse(ExAC_AMR < 0.01 | is.na(ExAC_AMR), TRUE, FALSE),
      TRUE ~ FALSE # Default case if none of the above conditions are met
    )) %>%
    filter(ExAC_Filtered) 
  
  # Further filtering based on ExonicFunc and Func
  inner1 <- binded1 %>% filter(
    !str_detect(ExonicFunc, "nonframeshift") | 
      str_detect(Func, "splicing")
  )
  
  inner2 <- inner1 %>% filter(
    str_detect(ExonicFunc, "frameshift") | 
      str_detect(ExonicFunc, "stop") | 
      str_detect(Func, "splicing")
  )
  
  inner3 <- binded1 %>% filter(
    REVEL >= 0.7 | str_detect(CLNDSDBID, "Pathogenic") | 
      str_detect(CLNDSDBID, "Likely")
  )
  
  # Combine inner2 and inner3
  inner_all <- rbind(inner2, inner3)
  
  # Filter out benign annotations
  inner_all <- inner_all %>% filter(
    !str_detect(CLNDSDBID, "benign") | is.na(CLNDSDBID)
  )
  
  inner_all <- inner_all %>% filter(
    !str_detect(CLNDSDBID, "Benign") | is.na(CLNDSDBID)
  )
  
  # Append to raw_variant
  raw_variant <- rbind(raw_variant, inner_all)
  
  # Print the current iteration
  print(i)
}

raw_variant <-raw_variant %>% distinct(germline_id,Gene,Exon,.keep_all = TRUE) #14713  # same gene, same exon -> remove
raw_variant3<-raw_variant[,c("germline_id","Gene","ExonicFunc","GenbankID","Nucleotide")] 
raw_variant_tmp<-dplyr::summarise(group_by(raw_variant3,Gene,ExonicFunc,GenbankID,Nucleotide),length(germline_id))


colnames(raw_variant_tmp)[which(names(raw_variant_tmp) == "length(germline_id)")] <- "counts_dup"

raw_variant2<-left_join(raw_variant,raw_variant_tmp, by=c("Gene","ExonicFunc","GenbankID","Nucleotide"))


raw_variant2<-raw_variant2 %>% relocate(counts_dup) #14713
raw_variant2$counts_dup<-as.numeric(raw_variant2$counts_dup)
raw_variant2<-raw_variant2 %>% filter(counts_dup<22) #9932

raw_variant3<-raw_variant2
nrow(distinct(raw_variant3,germline_id)) #

distinct_germline_count <- nrow(distinct(raw_variant3, germline_id))
gvb<- as.data.frame(table(raw_variant3$germline_id))
summary_stats <- summary(gvb$Freq) 
distinct_germline_count
summary_stats
summary(TARGET_data$germ_rare_burden)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.00   34.00   39.50   45.15   48.00  190.00 

colnames(gvb) <-c('germline_id','germ_rare_burden')
raw_variant3 <- mutate_at(raw_variant3, "ExonicFunc", ~replace(., is.na(.), "splicing"))


# 5. TMB --------------------------------------------------------------------------------------

somatic_mafs=list.files(path =paste0(data_folder,"/MAF/", pattern = "merged", full.names = TRUE))

merge_somatic=merge_mafs(mafs=somatic_mafs)

variant_per_sample<-merge_somatic@variants.per.sample

merge_somatic_qc<-merge_somatic
merge_data<-merge_somatic_qc@data #205294
merge_data<-merge_data %>% filter(t_depth>=30 & t_alt_count>2) #196602
merge_data<-merge_data %>% filter(t_alt_count/t_depth>0.05) #86657

merge_data2<-merge_data %>%
  filter(ifelse(Reference_Allele=="G"&Tumor_Seq_Allele2 %in% c("T"), t_alt_count/t_depth>=0.15, t_alt_count/t_depth>0.05)) %>%
  filter(ifelse(Reference_Allele=="C"&Tumor_Seq_Allele2 %in% c("A"), t_alt_count/t_depth>=0.15, t_alt_count/t_depth>0.05)) #11714


merge_data2<-merge_data2 %>% filter(is.na(gnomAD_AF)|gnomAD_AF<0.001) #10270


somatic_syn<-merge_somatic_qc@maf.silent #270110
somatic_syn<-somatic_syn %>% filter(t_depth>=30 & t_alt_count>2) #194406
somatic_syn<-somatic_syn %>% filter(t_alt_count/t_depth>0.05) #106248

somatic_syn<-somatic_syn %>% filter(is.na(gnomAD_AF)|gnomAD_AF<0.001) #97815

somatic_syn2<-somatic_syn %>%
  filter(ifelse(Reference_Allele=="G"&Tumor_Seq_Allele2 %in% c("T"), t_alt_count/t_depth>=0.15, t_alt_count/t_depth>0.05)) %>%
  filter(ifelse(Reference_Allele=="C"&Tumor_Seq_Allele2 %in% c("A"), t_alt_count/t_depth>=0.15, t_alt_count/t_depth>0.05)) #27834

colnames(merge_data2)[which(names(merge_data2) == "Tumor_Sample_Barcode")] <- "somatic_id"
colnames(somatic_syn2)[which(names(somatic_syn2) == "Tumor_Sample_Barcode")] <- "somatic_id"

all_data<-rbind(merge_data2,somatic_syn2)
Target_somatic_ID<-TARGET_data[,c("Tumor_Sample_Barcode","somatic_id")]
all_data<-left_join(all_data,Target_somatic_ID,by="somatic_id")
all_data<-all_data %>% relocate("Tumor_Sample_Barcode")
somatic.maf<-read.maf(maf =all_data)

TMB.data <- reshape2::melt(as.data.frame(getSampleSummary(somatic.maf))) %>%
  dplyr::mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels=getSampleSummary(somatic.maf)[[1]])) 

TMB_total<-TMB.data %>% filter(variable=="total")
TMB_total<-TMB_total[,c("Tumor_Sample_Barcode","value")]
colnames(TMB_total)<-c("Tumor_Sample_Barcode","TMB")

summary(TMB_total$TMB)
summary(TARGET_data$TMB_all)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   16.00   26.00   46.68   45.00  753.00 


# 6. Figure 2 ---------------------------------------------------------------------------------

colnames(TMB_total) <-c('Tumor_Sample_Barcode','TMB')

infor <- TARGET_data[,c("germline_id",'somatic_id','Tumor_Sample_Barcode')]

infor <-left_join(infor,TMB_total,by="Tumor_Sample_Barcode")
infor <-left_join(infor,gvb,by="germline_id")
infor$TMB_log <- log10(infor$TMB)

infor$germ_burden_group <- ifelse(infor$germ_rare_burden>=41,"High","Low")
cor.test(infor$germ_rare_burden,infor$TMB_log,method = 'spearman') # P = 0.018 : #Fig.1d
t.test(infor$TMB~infor$germ_burden_group) # 



TARGET_TMB_burden<-ggscatter(infor, x = "germ_rare_burden", y = "TMB_log",
                             add = "reg.line",                                 # Add regression line
                             conf.int = TRUE,
                             size=2,                                  # Add confidence interval
                             add.params = list(color = "black",
                                               fill = "lightgray"),
                             rug=FALSE,
                             fill= "#4AACC6",shape=21,color="#0E7AB0",
                             xlab="Germline variant burden of pFGVs",alpha=1
)+
  # stat_cor(method = "pearson",p.accuracy = 0.001, r.accuracy = 0.01,cor.coef.name='rho',
  #                    aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~")))+ 
  stat_cor(method = "spearman",
           p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.5, label.y.npc = 0.95,
           cor.coef.name='rho',
           aes(label = paste(
             ..r.label.., gsub("p", "P", ..p.label..), 
             sep = "~`,`~")))+ sm_hgrid()+
  ylab(expression("Somatic mutation burden (log"[10]*" scale)"))
TARGET_TMB_burden




# 7. DDR --------------------------------------------------------------------------------------
target_variant<-raw_variant3
ddr_variant<-target_variant %>% filter(Gene %in% DDR$Gene)
ddr_burden<-table(ddr_variant$germline_id)
ddr_burden<-data.frame(ddr_burden)
names(ddr_burden)<-c("germline_id","ddr_burden4")

TARGET_data <- left_join(TARGET_data,ddr_burden,by=c("germline_id"))
TARGET_data <- mutate_at(TARGET_data, "ddr_burden4", ~replace(., is.na(.), 0))
summary(TARGET_data$ddr_burden)
summary(TARGET_data$ddr_burden4)

TARGET_data <-  TARGET_data %>% mutate(DDR=ifelse(ddr_burden3>0,"Yes","No"))

TARGET_data2<-TARGET_data %>% filter(DDR=="No")
cor.test(TARGET_data2$TMB_log,TARGET_data2$germ_rare_burden,method = "spearman")

TMB_burden_target2<-ggscatter(TARGET_data2, x = "germ_rare_burden", y = "TMB_log",
                              add = "reg.line",                                 # Add regression line
                              conf.int = TRUE,
                              size=2,                                  # Add confidence interval
                              add.params = list(color = "black",
                                                fill = "lightgray"),
                              rug=FALSE,
                              fill= sm_color("skyblue"),shape=21,color=sm_color("blue"),stroke=0.5,
                              xlab="Germline variant burden of pFGVs",alpha=1
)+stat_cor(method = "spearman",p.accuracy = 0.001, r.accuracy = 0.01,cor.coef.name='rho',
           aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~")))+ sm_hgrid()+
  ylab(expression("Somatic mutation burden (log"[10]*" scale)"))



# 8. Z score ----------------------------------------------------------------------------------
TARGET_data$germ_z_score <- scale(TARGET_data$germ_rare_burden)
TARGET_data$outlier_excluded <- ifelse(TARGET_data$germ_z_score<3,"Normal",'Outlier')

df_filtered <- TARGET_data[abs(TARGET_data$germ_z_score) < 3, ]
cor.test(df_filtered$TMB_log,df_filtered$germ_rare_burden,method = "spearman")


# 9. Linear regression ------------------------------------------------------------------------

res <- lm(TARGET_data$TMB_log~TARGET_data$germ_rare_burden+TARGET_data$MEDIAN_Coverage_tumor+TARGET_data$MEDIAN_Coverage_germline+TARGET_data$Race)
summary(res)


# 10. Survival --------------------------------------------------------------------------------

## OS - only OS is reported
summary(TARGET_data$germ_rare_burden)
TARGET_data <- TARGET_data %>% mutate(germ_burden_group=ifelse(
  germ_rare_burden>=45,"High","Low"))
TARGET_data$germ_burden_group <- factor(TARGET_data$germ_burden_group,levels=c("Low","High"))

fit3<-survfit(Surv(fu_overall,event_death=="Dead")~germ_burden_group,data=TARGET_data)

ggsurvplot(fit3,data=TARGET_data,risk.table = TRUE,pval=TRUE)
sd<-survdiff(Surv(fu_overall, event_death=="Dead") ~ germ_burden_group, data = TARGET_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # 0.0051

## Multi Cox
res.cox <- coxph(Surv(fu_overall, event_death=="Dead") ~ germ_burden_group+MYCN+`Estimated population`, data=TARGET_data)
summary(res.cox)


# 11. CPG -------------------------------------------------------------------------------------

target_pfgv <- read_tsv(paste0(output_folder,'/TARGET_gHFV_All.tsv'))


setdiff(raw_variant3$Gene,raw_variant4$Gene)


## Cosmic All
cosmic<-read_csv(paste0(data_folder,'/ref_gene/cancer_gene_census.csv'))

cosmic_TSG_onco<- cosmic %>% filter(str_detect(Role,'TSG|suppresses|oncogene'))
cosmic_CPG<-cosmic_TSG_onco %>% filter(Germline=="yes")
cosmic_TSG<-cosmic %>% filter(str_detect(Role,'TSG|suppresses'))
cosmic_TSG_Only<-cosmic_TSG %>% filter(!str_detect(Role,'onco'))

cosmic_TSG_germline<-cosmic_TSG %>%  filter(Germline=="yes")


missense<-target_pfgv %>% filter(ExonicFunc=="nonsynonymous SNV") %>% filter(Gene %in% cosmic_TSG_onco$Gene) #107
PTV<-target_pfgv %>% filter(ExonicFunc!="nonsynonymous SNV") %>% filter(Gene %in% cosmic_TSG_Only$Gene) #121

Cosmic<-rbind(missense,PTV)

Cosmic<-Cosmic %>% filter(germline_id %in% TARGET_data$germline_id)

cosmic_burden<-table(Cosmic$germline_id)
cosmic_burden<-data.frame(cosmic_burden)
names(cosmic_burden)<-c("germline_id","cosmic_burden")

TARGET_data <- left_join(TARGET_data,cosmic_burden,by=c("germline_id"))
TARGET_data <- mutate_at(TARGET_data, "cosmic_burden", ~replace(., is.na(.), 0))

Cosmic_CPG<-Cosmic %>% filter(Gene %in% cosmic_CPG$Gene)
nrow(Cosmic_CPG) #97
nrow(distinct(Cosmic_CPG,germline_id)) #79
table(TARGET_data$CPG) #79
length(unique(Cosmic_CPG$Gene)) #43


# CPG Survival --------------------------------------------------------------------------------
TARGET_data_non <- TARGET_data %>% filter(mycn2!="Amplified")
TARGET_data_amp <- TARGET_data %>% filter(mycn2=="Amplified")

fit<-survfit(Surv(fu_overall,event_death=="Dead")~ CPG,data=TARGET_data)


sd<-survdiff(Surv(fu_overall, event_death=="Dead") ~ CPG, data = TARGET_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2

fit2<-survfit(Surv(fu_overall,event_death=="Dead")~ CPG,data=TARGET_data_non)
fit2

sd<-survdiff(Surv(fu_overall, event_death=="Dead") ~ CPG, data = TARGET_data_non)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2

fit3<-survfit(Surv(fu_overall,event_death=="Dead")~ CPG,data=TARGET_data_amp)

sd<-survdiff(Surv(fu_overall, event_death=="Dead") ~ CPG, data = TARGET_data_amp)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2


