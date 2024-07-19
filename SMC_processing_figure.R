
# 1. Library ----------------------------------------------------------------------------------
package_list <- c("readr", 'tidyverse',"tibble", "openxlsx", "readxl", 
                  "writexl", "stringr", "naniar", "MASS", "devtools", 
                  "ggvenn", "cowplot", "ggpubr", "survival", "survminer",
                  "hrbrthemes", "viridis", "RColorBrewer", "scales", "psych",'smplot2')
purrr::walk(package_list, library, character.only = TRUE)

output_folder <- '/path/output/SMC'
data_folder <-'/path/data/SMC/'
# 2. Data -------------------------------------------------------------------------------------

SMC_data <- readRDS(paste0(data_folder,'/SMC_data.Rds'))

# 3. Germline Functional ------------------------------------------------------------------

allfilelist <- list.files(path = data_folder, pattern = "smc", recursive = FALSE, full.names = TRUE)


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
  tmpdata$KRGDB_1100 <- as.numeric(tmpdata$KRGDB_1100)
  tmpdata$KRGDB_622 <- as.numeric(tmpdata$KRGDB_622)
  
  # Extract resource ID info from the filename
  resource_id_info <- str_sub(allfilelist[i], 62, 78)
  resource_id_column <- data.frame(germline_id = resource_id_info)
  
  # Apply filters to the data
  tmpdata <- tmpdata %>% 
    filter(DP >= 10) %>% 
    filter(GQ >= 50) %>% 
    filter(VAF >= 0.3) %>% 
    filter(!str_detect(FILTER, "PASS")) %>%  #VQSRTrancheSNP99.9
    filter(FILTER != ".") %>% 
    filter(Chr != "chrY")
  
  # Combine resource ID with the filtered data
  binded1 <- cbind(resource_id_column, tmpdata)
  
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

raw_variant <- mutate_at(raw_variant, "ExonicFunc", ~replace(., is.na(.), "splicing"))

# Remove duplicates
raw_variant <- raw_variant %>% distinct(germline_id, Gene, Exon, .keep_all = TRUE) # same gene, same exon -> remove (duplication)

# Select relevant columns and count duplicates
raw_variant3 <- raw_variant[, c("germline_id", "Gene", "ExonicFunc", "GenbankID", "Nucleotide")]
raw_variant_tmp <- raw_variant3 %>% group_by(Gene, ExonicFunc, GenbankID, Nucleotide) %>% dplyr::summarise(counts_dup = n())
raw_variant2 <- left_join(raw_variant, raw_variant_tmp, by = c("Gene", "ExonicFunc", "GenbankID", "Nucleotide"))

# Relocate counts_dup column
raw_variant2 <- raw_variant2 %>% relocate(counts_dup)

# Filter based on counts_dup
raw_variant3 <- raw_variant2 %>% filter(counts_dup < 13)

# write_tsv(raw_variant3,paste0(output_folder,'/SMC_gHFV_All.tsv'))

# Count distinct germline_id and create a summary
distinct_germline_count <- nrow(distinct(raw_variant3, germline_id))
gvb<- as.data.frame(table(raw_variant3$germline_id))
summary_stats <- summary(gvb$Freq) 
distinct_germline_count
summary_stats
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 27.00   37.00   41.00   41.58   46.00   58.00 

colnames(gvb) <-c('germline_id','germ_rare_burden')

ptv_all<-raw_variant3 %>% filter(ExonicFunc!="nonsynonymous SNV"|is.na(ExonicFunc))

ptv_burden<-table(ptv_all$germline_id)
ptv_burden<-data.frame(ptv_burden)
names(ptv_burden)<-c("Tumor_Sample_Barcode","germ_ptv_burden")
summary(ptv_burden$germ_ptv_burden)
summary(SMC_data$germ_ptv_burden)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13.00   22.00   24.00   25.08   28.00   39.00 
# write_tsv(raw_variant3,paste0(output_folder,'/SMC_gHFV_All.tsv'))

raw_variant3 <-read_tsv(paste0(output_folder,'/SMC_gHFV_All.tsv'))

# 4. Somatic Mutation ------------------------------------------------------------------

somatic_mafs=list.files(path =paste0(data_folder,"/MAF/", pattern = "merged", full.names = TRUE))
merge_somatic=merge_mafs(mafs=somatic_mafs)

variant_per_sample<-merge_somatic@variants.per.sample

merge_somatic_qc<-merge_somatic
merge_data<-merge_somatic_qc@data 
merge_data<-merge_data %>% filter(t_depth>=30 & t_alt_count>2)
merge_data<-merge_data %>% filter(t_alt_count/t_depth>0.05)
merge_data<-merge_data %>% filter(is.na(gnomAD_EAS_AF)|gnomAD_EAS_AF<0.001)

somatic_syn<-merge_somatic_qc@maf.silent 
somatic_syn<-somatic_syn %>% filter(t_depth>=30 & t_alt_count>2) #
somatic_syn<-somatic_syn %>% filter(t_alt_count/t_depth>0.05) #
somatic_syn<-somatic_syn %>% filter(is.na(gnomAD_EAS_AF)|gnomAD_EAS_AF<=0.001)

colnames(merge_data)[which(names(merge_data) == "Tumor_Sample_Barcode")] <- "somatic_id"
colnames(somatic_syn)[which(names(somatic_syn) == "Tumor_Sample_Barcode")] <- "somatic_id"

all_data<-rbind(merge_data,somatic_syn)
SMC_somatic_ID<-SMC_data[,c("Tumor_Sample_Barcode","somatic_id")]
all_data<-left_join(all_data,SMC_somatic_ID,by="somatic_id")
all_data<-all_data %>% relocate("Tumor_Sample_Barcode")
somatic.maf<-read.maf(maf =all_data)

somatic_all.maf=subsetMaf(maf=somatic.maf,
                          clinQuery="Tumor_Sample_Barcode %in% SMC_data$Tumor_Sample_Barcode")

TMB.data <- reshape2::melt(as.data.frame(getSampleSummary(somatic.maf))) %>%
  # dplyr::inner_join(TARGET_clinic) %>%
  dplyr::mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels=getSampleSummary(somatic.maf)[[1]])) 

# TARGET_clinic_white <- TARGET_clinic %>% filter(Race=="White")
TMB_total<-TMB.data %>% filter(variable=="total")
TMB_total<-TMB_total[,c("Tumor_Sample_Barcode","value")]
colnames(TMB_total)<-c("Tumor_Sample_Barcode","TMB")
summary(TMB_total$TMB)
summary(SMC_data$TMB)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00   10.00   17.00   23.74   30.00  143.00 

# write_tsv(TMB.data,paste0(output_folder,'/SMC_TMB_data.tsv'))

# 5. Germline Syn ------------------------------------------------------------------


allfilelist = list.files(path =data_folder, pattern = "synonymous", recursive = TRUE, full.names = TRUE)
filecnt = length(allfilelist)
germ_syn <- data.frame()
for (i in 1:filecnt) {
  tmpdata <- read.delim(allfilelist[i])
  
  tmpdata$n_count_alt<-as.numeric(tmpdata$n_count_alt)
  tmpdata$DP<-as.numeric(tmpdata$DP)
  tmpdata<-tmpdata %>% mutate(VAF=n_count_alt/DP)
  tmpdata$VAF<-as.numeric(tmpdata$VAF)
  tmpdata$GQ<-as.numeric(tmpdata$GQ)
  tmpdata$AF<-as.numeric(tmpdata$AF)
  tmpdata<-tmpdata %>% filter(VAF>=0.3)
  tmpdata<-tmpdata %>% filter(GQ>=50)
  tmpdata<-tmpdata %>% filter(DP>=10)
  tmpdata<-tmpdata %>% filter(KRGDB_1100<0.01|is.na(KRGDB_1100))
  tmpdata<-tmpdata %>% filter(ExAC_ALL<0.01|is.na(ExAC_ALL))
  tmpdata<-tmpdata %>% filter(!str_detect(FILTER,"100")) 
  tmpdata<-tmpdata %>% filter(Chr!="chrY") 
  tmpdata$KRGDB_1100<-as.numeric(tmpdata$KRGDB_1100)
  tmpdata$KRGDB_622<-as.numeric(tmpdata$KRGDB_622)
  tmpdata$ExAC_ALL<-as.numeric(tmpdata$ExAC_ALL) 
  resource_id_info <- str_sub(allfilelist[i],77,93)
  resource_id_column<-data.frame(germline_id=resource_id_info)
  binded1 <- cbind(resource_id_column,tmpdata)
  germ_syn <-rbind(germ_syn,binded1)
  print(i)
}


germ_syn <- germ_syn %>% filter(germline_id %in% SMC_data$germline_id) #30336
germ_syn3<-germ_syn[,c("germline_id","Gene.refGene","Start","End","Ref","Alt")]

germ_syn3_tmp<-dplyr::summarise(group_by(germ_syn3,Gene.refGene,,Start,End,Ref,Alt),length(germline_id))
colnames(germ_syn3_tmp)[which(names(germ_syn3_tmp) == "length(germline_id)")] <- "counts_dup"

germ_syn3_tmp2<-dplyr::left_join(germ_syn3,germ_syn3_tmp, by=c("Gene.refGene","Start","End","Ref","Alt"))


germ_syn3_tmp2<-germ_syn3_tmp2 %>% relocate(counts_dup)
germ_syn3_tmp2$counts_dup<-as.numeric(germ_syn3_tmp2$counts_dup)
germ_syn<-left_join(germ_syn,germ_syn3_tmp2,by=c("germline_id","Gene.refGene","Start","End","Ref","Alt"))
germ_syn_filtered<-germ_syn %>% filter(counts_dup<13) #27280 
germ_syn_table <- as.data.frame(table(germ_syn_filtered$germline_id))
summary(germ_syn_table$Freq)
summary(SMC_data$germ_syn_burden)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 190.0   208.0   219.0   218.2   228.0   254.0 

# write_tsv(germ_syn_filtered,paste0(output_folder,'/SMC_Syn_Variant.tsv'))

# 6. DDR -----------------------------------------------------------------------------------------
DDR<-read_xlsx(paste0(data_folder,"/DDR_core.xlsx"))

smc_pfgv_variants <-read_tsv(paste0(output_folder,'/SMC_gHFV_All.tsv')) #5198

ddr_variant<-smc_pfgv_variants %>% filter(Gene %in% DDR$Gene)
ddr_burden<-table(ddr_variant$germline_id)
ddr_burden<-data.frame(ddr_burden)
names(ddr_burden)<-c("germline_id","ddr_burden")
summary(ddr_burden$ddr_burden)

SMC_data <- left_join(SMC_data,ddr_burden,by=c("germline_id"))
SMC_data <- mutate_at(SMC_data, "ddr_burden", ~replace(., is.na(.), 0))
SMC_data <-  SMC_data %>% mutate(DDR=ifelse(ddr_burden>0,"Yes","No"))
shapiro.test(SMC_data$TMB)
shapiro.test(SMC_data$TMB_log)
shapiro.test(SMC_data$germ_rare_burden)

SMC_data2<-SMC_data %>% filter(DDR=="No")
cor.test(SMC_data2$TMB_log,SMC_data2$germ_rare_burden,method = "pearson") #0.23, 0.032



# 7. Figure 1 --------------------------------------------------------------------------------

colnames(TMB_total) <-c('Tumor_Sample_Barcode','TMB')

infor <- SMC_data[,c("germline_id",'somatic_id','Tumor_Sample_Barcode')]

infor <-left_join(infor,TMB_total,by="Tumor_Sample_Barcode")
infor <-left_join(infor,gvb,by="germline_id")
infor$TMB_log <- log10(infor$TMB)
infor$germ_burden_group <- ifelse(infor$germ_rare_burden>=41,"High","Low")
cor.test(infor$germ_rare_burden,infor$TMB_log) # P = 0.018 : #Fig.1d
t.test(infor$TMB~infor$germ_burden_group) # 0.018

TMB_burden<-ggscatter(SMC_data, x = "germ_rare_burden", y = "TMB_log",
                      add = "reg.line",                                 # Add regression line
                      conf.int = TRUE,
                      size=2,                                  # Add confidence interval
                      add.params = list(color = "black",
                                        fill = "lightgray"),
                      rug=FALSE,
                      fill= sm_color("skyblue"),shape=21,color=sm_color("blue"),stroke=0.5,
                      xlab="Germline variant burden of pFGVs",alpha=1
)+stat_cor(method = "pearson",p.accuracy = 0.001, r.accuracy = 0.01,cor.coef.name='r',
           aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~")))+ sm_hgrid()+
  ylab(expression("Somatic mutational burden (log"[10]*" scale)"))
# ylab(expression("log"[10]*" Somatic mutational burden"))
TMB_burden

cor.test(SMC_data$TMB_log,SMC_data$germ_syn_burden,method="spearman")
TMB_germ_syn<-ggscatter(SMC_data, x = "germ_syn_burden", y = "TMB_log",
                        add = "reg.line",                                 # Add regression line
                        conf.int = TRUE,
                        size=2,                                  # Add confidence interval
                        add.params = list(color = "black",
                                          fill = "lightgray"),
                        rug=FALSE,shape=21,
                        fill= sm_color("skyblue"),color=sm_color("blue"),stroke=0.5,
                        xlab="Number of rare synonymous germline variant"
)+stat_cor(method = "pearson",p.accuracy = 0.001, r.accuracy = 0.01,cor.coef.name='r',
           aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~")))+ 
  ylab(expression("Somatic mutational burden (log"[10]*" scale)"))+sm_hgrid()

TMB_germ_syn_histogram<-ggExtra::ggMarginal(TMB_germ_syn, type = "histogram",fill="#CDCDCD")

TMB_burden2<-ggscatter(SMC_data2, x = "germ_rare_burden", y = "TMB_log",
                       add = "reg.line",                                 # Add regression line
                       conf.int = TRUE,
                       size=2,                                  # Add confidence interval
                       add.params = list(color = "black",
                                         fill = "lightgray"),
                       rug=FALSE,
                       fill= sm_color("skyblue"),shape=21,color=sm_color("blue"),stroke=0.5,
                       xlab="Germline variant burden of pFGVs",alpha=1
)+stat_cor(method = "pearson",p.accuracy = 0.001, r.accuracy = 0.01,cor.coef.name='r',
           aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~")))+ sm_hgrid()+
  ylab(expression("Somatic mutational burden (log"[10]*" scale)"))


TMB_burden_histogram2<-ggExtra::ggMarginal(TMB_burden2, type = "histogram",fill="#CDCDCD")


GVB_TMB_Box<-SMC_data %>% 
  ggplot(aes(x=germ_burden_group, fill= germ_burden_group, y=TMB))+
  geom_point(position=position_jitterdodge(),
             color="black",
             shape=21,
             # fill="gray",
             stroke=0.5)+
  geom_boxplot(show.legend = FALSE,outlier.shape = NA,alpha=c(0.7,0.7),width=0.5) +
  scale_fill_manual(values = sm_color("brown","blue"))+
  ggpubr::stat_compare_means(method = "t.test",
                             label.y=75,comparisons = my_comparisons) +
  xlab("Germline variant burden of pFGVs")+
  ylab("Somatic mutational burden")+theme(
    panel.background = element_blank(),
    axis.title = element_text(face="bold"),
    axis.text.x = element_text(face="bold",color = "black"),
    axis.line = element_line(linewidth =0.5, colour = "black"),
    legend.position="none")+sm_hgrid()+
  scale_x_discrete(labels=c('Low', 'High'))


# 8. Regression -------------------------------------------------------------------------------
res <- lm(SMC_data$TMB_log~SMC_data$germ_rare_burden+SMC_data$Germline_Coverage+SMC_data$Tumor_Coverage+SMC_data$Risk)
# SMC_data$germ_rare_burden   0.009339   0.004618   2.022  0.04537 *  


# 9. Survival ---------------------------------------------------------------------------------

summary(SMC_data$germ_rare_burden) #Mean 41 & Median also 41 
SMC_data <- SMC_data %>% mutate(germ_burden_group=ifelse(germ_rare_burden>=41,"High","Low"))

SMC_data$germ_burden_group<-factor(SMC_data$germ_burden_group,levels=c("Low","High"))

## PFS
fit<-survfit(Surv(fu_event_pr,event_pr=="Yes")~germ_burden_group,data=SMC_data)

ggsurvplot(fit,data=SMC_data,risk.table = TRUE,pval=TRUE)
sd<-survdiff(Surv(fu_event_pr, event_pr=="Yes") ~ germ_burden_group, data = SMC_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # 0.018 

## OS
fit2<-survfit(Surv(fu_overall,event_death=="Yes")~germ_burden_group,data=SMC_data)
fit2
ggsurvplot(fit2,data=SMC_data,risk.table = TRUE,pval=TRUE)
sd<-survdiff(Surv(fu_overall, event_death=="Yes") ~ germ_burden_group, data = SMC_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # 0.455

## Multi Cox
res.cox <- coxph(Surv(fu_event_pr, event_pr=="Yes") ~ germ_burden_group+Stage+mycn+Age2, data=SMC_data)
summary(res.cox)
# exp(coef) 2.781, P = 0.01838


# 10. CPG -------------------------------------------------------------------------------------
smc_pfgv <- read_tsv(paste0(output_folder,'/SMC_gHFV_All.tsv'))
cosmic<-read_csv(paste0(data_folder,'/ref_gene/cancer_gene_census.csv'))


## Cosmic All

cosmic_TSG_onco<- cosmic %>% filter(str_detect(Role,'TSG|suppresses|oncogene'))
cosmic_CPG<-cosmic_TSG_onco %>% filter(Germline=="yes")
cosmic_TSG<-cosmic %>% filter(str_detect(Role,'TSG|suppresses'))
cosmic_TSG_Only<-cosmic_TSG %>% filter(!str_detect(Role,'onco'))

cosmic_TSG_germline<-cosmic_TSG %>%  filter(Germline=="yes")

missense<-smc_pfgv %>% filter(ExonicFunc=="nonsynonymous SNV") %>% filter(Gene %in% cosmic_TSG_onco$Gene) #61
PTV<-smc_pfgv %>% filter(ExonicFunc!="nonsynonymous SNV") %>% filter(Gene %in% cosmic_TSG_Only$Gene) #33


Cosmic<-rbind(missense,PTV)
Cosmic<-Cosmic %>% filter(germline_id %in% SMC_data$germline_id)

cosmic_burden<-table(Cosmic$germline_id)
cosmic_burden<-data.frame(cosmic_burden)
names(cosmic_burden)<-c("germline_id","cosmic_burden")

SMC_data <- left_join(SMC_data,cosmic_burden,by=c("germline_id"))
SMC_data <- mutate_at(SMC_data, "cosmic_burden", ~replace(., is.na(.), 0))

nrow(distinct(Cosmic,Gene)) # 70 
nrow(distinct(Cosmic,germline_id)) #68

summary(SMC_data$cosmic_burden)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   1.000   0.752   1.000   3.000 

## CPG
Cosmic_CPG<-Cosmic %>% filter(Gene %in% cosmic_CPG$Gene)
nrow(Cosmic_CPG) #45
nrow(distinct(Cosmic_CPG,germline_id)) #39
table(SMC_data$CPG) #39 
length(unique(Cosmic_CPG$Gene)) #31


## CPG clinical Factor -------------------------------------------------------------------------
library(gtsummary)
library(webshot2)


SMC_data %>% 
  dplyr::select(CPG,Age2,age, sex,fam_2nd_bi,stage,sex,path_class, mycn,site, risk, Risk,event_pr, event_death, TMB, germ_rare_burden) %>% 
  tbl_summary(by = CPG, type = all_dichotomous() ~ "categorical",
              statistic = list(all_continuous() ~ "{mean}",
                               all_categorical() ~ "{n} / {N} ({p}%)"),
              digits = all_continuous() ~ 2) %>% 
  add_n() %>% # add column with total number of non-missing observations
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 3)) %>% # test for a difference between groups
  modify_header(label = "**Characteristics**") %>% # update the column header
  bold_labels() 


## ACMG ----------------------------------------------------------------------------------------

ACMG<-read_delim(paste0(data_folder,'/SMC_ACMG_Annot.txt'))
ACMG<-ACMG[,c("germline_id","Gene","ExonicFunc","ACMG evidence","GenbankID","AA","ACMG")]

CPG_ACMG<-left_join(Cosmic_CPG,ACMG,by=c("germline_id","Gene","ExonicFunc","GenbankID","AA"))
CPG_ACMG


ACMG_table<-table(CPG_ACMG$germline_id,CPG_ACMG$ACMG)

ACMG_data<-data.frame(rbind(ACMG_table))
ACMG_data$germline_id<-rownames(ACMG_data)
rownames(ACMG_data)<-c(1:nrow(ACMG_data))
ACMG_data<-ACMG_data %>% mutate(PV_LPV=ifelse(PV==1|LPV==1,"Yes","No"))
ACMG_data_PV<-ACMG_data %>% filter(PV_LPV=="Yes")

acmg <- unique(ACMG_data_PV$germline_id)  
length(acmg) #16



ACMG_data_VUS<-ACMG_data %>% filter(PV_LPV=="No")

SMC_data<-SMC_data %>% mutate(CPG=ifelse(germline_id %in% CPG$germline_id,"Yes","No"))
SMC_data<-SMC_data %>% mutate(ACMG=ifelse(germline_id %in% ACMG_data_PV$germline_id,"PV/LPV",
                                              ifelse(germline_id %in% ACMG_data_VUS$germline_id,"VUS","No")))

SMC_data<-SMC_data %>% mutate(ACMG2=ifelse(ACMG=="PV/LPV","PV/LPV","No"))


## CPG Survival --------------------------------------------------------------------------------

## PFS 
fit<-survfit(Surv(fu_event_pr,event_pr=="Yes")~CPG,data=SMC_data)
fit
ggsurvplot(fit,data=SMC_data,risk.table = TRUE,pval=TRUE)

sd<-survdiff(Surv(fu_event_pr, event_pr=="Yes") ~ CPG, data = SMC_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # <0.001

## OS 
fit<-survfit(Surv(fu_overall,event_death=="Yes")~CPG,data=SMC_data)
fit
ggsurvplot(fit,data=SMC_data,risk.table = TRUE,pval=TRUE)

sd<-survdiff(Surv(fu_overall, event_death=="Yes") ~ CPG, data = SMC_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # 0.0249


## ACMG Survival --------------------------------------------------------------------------------

## PFS 
fit<-survfit(Surv(fu_event_pr,event_pr=="Yes")~ACMG2,data=SMC_data)
fit
ggsurvplot(fit,data=SMC_data,risk.table = TRUE,pval=TRUE)

sd<-survdiff(Surv(fu_event_pr, event_pr=="Yes") ~ ACMG2, data = SMC_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # <0.001

## OS 
fit<-survfit(Surv(fu_overall,event_death=="Yes")~ACMG2,data=SMC_data)
fit
ggsurvplot(fit,data=SMC_data,risk.table = TRUE,pval=TRUE)

sd<-survdiff(Surv(fu_overall, event_death=="Yes") ~ ACMG2, data = SMC_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # 0.009



## CPG Uni Cox  ------------------------------------------------------------------------------------

library(moonBook)
library(forestmodel)

SMC_data_tmp<-SMC_data[,c("Stage","Age2","CPG","gHFV","germ_burden_group","mycn","Risk","fu_event_pr","event_pr","fu_overall","event_death")]
SMC_data_tmp$mycn<-factor(SMC_data_tmp$mycn,levels=c("non-MNA","MNA"))
SMC_data_tmp$MYCN<-factor(SMC_data_tmp$mycn,levels=c("non-MNA","MNA"),labels=c("Not amplified","Amplified"))
SMC_data_tmp$Age<-factor(SMC_data_tmp$Age2)
SMC_data_tmp$Risk<-factor(SMC_data_tmp$Risk,levels=c("non-high","high"),labels=c("Non-high risk","High risk"))

SMC_data_tmp$germ_burden_group<-factor(SMC_data_tmp$germ_burden_group,levels=c("Low","High"),labels = c("Low","High"))
SMC_data_tmp$`Germline variant burden`<-SMC_data_tmp$germ_burden_group
SMC_data_tmp$`pFGVs in CPGs`<-SMC_data_tmp$CPG
SMC_data_tmp$`pFGVs in cancer-relevant genes`<-SMC_data_tmp$gHFV

vars_for_table<-c("Age","Stage","MYCN","Risk","`Germline variant burden`","`pFGVs in CPGs`")
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(fu_event_pr, event_pr=="Yes")~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = SMC_data_tmp)})


A<-forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)

pdf(paste0(output_folder, "/SMC_cox_univariate_empty.pdf"), width=10, height=4,onefile = FALSE)
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T,breaks=c(0,log(2),log(5),log(10),log(30)),
                   ,limits=c(-log(2),log(50)),panels = panels))
dev.off()

pdf(paste0(output_folder, "/SMC_cox_univariate.pdf"), width=10, height=4,onefile = FALSE)
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T,breaks=c(0,log(2),log(5),log(10),log(30)),
                   ,limits=c(-log(2),log(50))))
dev.off()

panels <- list( list(width = 0.03), list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"), list(width = 0.1, display = ~level), list(width = 0.05, display = ~n, hjust = 1, heading = "N"), 
                list(width = 0.03, item = "vline", hjust = 0.5), list( width = 0.55, item = "forest", hjust = 0.5, heading = "Hazard ratio", linetype = "dashed", line_x = 0 ), list(width = 0.03, item = "vline", hjust = 0.5), list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf( "%0.2f (%0.2f, %0.2f)", 2, 2, 2)), display_na = NA), list( width = 0.05, display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.1)), display_na = NA, hjust = 1, heading = "p" ), list(width = 0.03) )


## Cox Multi -----------------------------------------------------------------------------------

### PFS
SMC_multi <- SMC_data_tmp %>%
  transmute(fu_event_pr,
            event_pr,
            Age = Age,
            Stage = Stage,
            MYCN=MYCN,
            Risk=Risk,
            `Germline variant burden`=`Germline variant burden`,
            `pFGVs in CPGs` = `pFGVs in CPGs`
  )


print(forest_model(coxph(Surv(fu_event_pr, event_pr=="Yes") ~ ., SMC_multi),breaks=c(0,log(2),log(5),log(10),log(30)),
                   limits=c(-log(2),log(50))))

res.cox <- coxph(Surv(fu_event_pr, event_pr=="Yes") ~ Age+Stage+MYCN+Risk+`Germline variant burden`+`pFGVs in CPGs`, data=SMC_data_tmp)
summary(res.cox)
#`pFGVs in CPGs`Yes     HR 2.9089  P = 0.0065 **


SMC_data_nocpg<-SMC_data %>% filter(CPG=="No")
sd<-survdiff(Surv(fu_event_pr, event_pr=="Yes") ~ germ_burden_group, data = SMC_data_nocpg)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 #0.0127



# MYCN ----------------------------------------------------------------------------------------
table(SMC_data$mycn)
nonMNA<-SMC_data %>% filter(mycn=="non-MNA")
MNA<-SMC_data %>% filter(mycn=="MNA")

table(SMC_data$mycn)

fit_os<-survfit(Surv(fu_overall,event_death=="Yes")~CPG,data=nonMNA)

ggsurvplot(fit_os,data=nonMNA,risk.table = TRUE,pval=TRUE)

sd<-survdiff(Surv(fu_overall, event_death=="Yes") ~ CPG, data = nonMNA)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # 0.050


sd<-survdiff(Surv(fu_overall, event_death=="Yes") ~ CPG, data = MNA)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2 # 0.223


# Cosmic --------------------------------------------------------------------------------------

SMC_data<-SMC_data %>% mutate(Cosmic=ifelse(CPG=="Yes","CPG",
                                                ifelse(gHFV=="Yes","Others","No")))
table(SMC_data$Cosmic)
SMC_data$Cosmic<-factor(SMC_data$Cosmic,levels=c("No","Others","CPG"))

sd<-survdiff(Surv(fu_event_pr, event_pr=="Yes") ~ Cosmic, data =SMC_data)
p_val<-1 - pchisq(sd$chisq, length(sd$n) - 1)
p_val2<-round(p_val,4)
p_val2

