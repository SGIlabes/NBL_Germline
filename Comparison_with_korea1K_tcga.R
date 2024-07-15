
# 1. Library ----------------------------------------------------------------------------------
package_list <- c("readr", 'tidyverse',"tibble", "openxlsx", "readxl", 
                  "writexl", "stringr", "naniar", "MASS", "devtools", 
                  "ggvenn", "cowplot", "ggpubr", "survival", "survminer",
                  "hrbrthemes", "viridis", "RColorBrewer", "scales", "psych",'smplot2')
purrr::walk(package_list, library, character.only = TRUE)

output_folder <- '/path/output/'
data_folder <-'/path/data//'

TARGET_data <- readRDS(paste0(data_folder,'/TARGET_data.RDS'))
SMC_data <- readRDS(paste0(data_folder,'/SMC_data.RDS'))
TCGA_data <- readRDS(paste0(data_folder,'/TCGA_data.RDS'))
KOREA1K <- readRDS(paste0(data_folder,'/KOREA1K_data.RDS'))


# 2. Enrichment -------------------------------------------------------------------------------------

## TARGET vs. TCGA
TARGET_data_white <- TARGET_data %>% filter(Race=="White")

TARGET_tmp <- TARGET_data_white[,c("Tumor_Sample_Barcode",'CPG')]
TARGET_tmp$cat <- 'TARGET'
TCGA_tmp <- TCGA_data[,c("Tumor_Sample_Barcode",'CPG')]
TCGA_tmp$cat <-'TCGA'

all <- rbind(TARGET_tmp,TCGA_tmp)
all$cat <- factor(all$cat,levels=c("TCGA",'TARGET'))
fisher.test(all$cat,all$CPG)

## KOREA1K vs. SMC


SMC_data_tmp <- SMC_data[,c("Tumor_Sample_Barcode",'CPG')]
SMC_data_tmp$cat <- 'SMC'
KOREA1K_tmp <- KOREA1K[,c("Tumor_Sample_Barcode",'CPG')]
KOREA1K_tmp$cat <-'Korea1K'

all <- rbind(SMC_data_tmp,KOREA1K_tmp)
all$cat <- factor(all$cat,levels=c("SMC",'Korea1K'))
fisher.test(all$cat,all$CPG)

# 4. Trend cor -----------------------------------------------------------------------------------

TCGA_data<-TCGA_data %>% mutate(Age_group=case_when(age_dx<30 ~ "Under 30",
                                                                                  age_dx<40 ~ "30s",
                                                                                  age_dx<50 ~ "40s",
                                                                                  age_dx<60 ~ "50s",
                                                                                  age_dx<70 ~ "60s",
                                                                                  age_dx<80 ~ "70s",
                                                                                  age_dx>=80 ~"Over 80"))

TCGA_tmp <- TCGA_data[,c("Tumor_Sample_Barcode","Age_group","germ_rare_burden","TMB")]
#TARGET_data <- readRDS(paste0(rds_dir,"/TARGET_data_Final.Rds"))
TARGET_white_tmp <- TARGET_data_white[,c("Tumor_Sample_Barcode","germ_rare_burden","TMB_all")]
colnames(TARGET_white_tmp)[which(names(TARGET_white_tmp) == "TMB_all")] <- "TMB"
TARGET_white_tmp$Age_group <-"TARGET-NBL"
TARGET_white_tmp$group <-"TARGET-NBL"
TCGA_tmp$group <-"TCGA"
All_age <- rbind(TCGA_tmp,TARGET_white_tmp)

All_age$Age_group<-factor(All_age$Age_group,levels=c("TARGET-NBL","Under 30","30s","40s","50s","60s","70s","Over 80"),
                          label=c("<18","19-29","30-39","40-49","50-59","60-69","70-79","Over 80"))

coeff_age<-All_age %>% group_by(Age_group) %>% group_modify(~ glance(spearman.ci(.x$germ_rare_burden, .x$TMB)))
coeff_age <- coeff_age %>% mutate(group=ifelse(Age_group=="<18","TARGET-NBL","TCGA"))
coeff_age$All<-"All"
annotate(geom = "text",
         x = 1:nrow(data),
         y = min(data$y),
         label = data$axis1,
         vjust = 3.5) +
  annotate(geom = "text",
           x = 1:nrow(data),
           y = min(data$y),
           label = data$axis2,
           vjust = 5)

library(grid)
PMCMRplus::jonckheereTest(coeff_age$estimate,coeff_age$Age_group,alternative = 'less')


# 5. CPG Trend -----------------------------------------------------------------------------------
library(janitor)
TARGET_tmp <- TARGET_data_white[,c("Tumor_Sample_Barcode","CPG")]
TARGET_tmp$Age_group <- "TARGET"
TCGA_tmp <- TCGA_data[,c("Tumor_Sample_Barcode","CPG","Age_group")]

All_white <- rbind(TARGET_tmp,TCGA_tmp)
All_white$Age_group<-factor(All_white$Age_group,levels=c("TARGET","Under 30","30s","40s","50s","60s","70s","Over 80"))

all_white_group<-All_white %>% group_by(Age_group,CPG) %>% dplyr::summarise(total_count=n()) %>%as.data.frame()

All_white$CPG <-factor(All_white$CPG, levels=c("Yes","No"))


count_table <- as.data.frame(All_white %>%
                               tabyl(Age_group, CPG) %>% adorn_totals(c("col")))

cpg_pts  <- count_table$Yes
all_pts <-count_table$Total
prop.test(cpg_pts, all_pts)
prop.trend.test(cpg_pts, all_pts)


freq_table<-All_white %>% 
  tabyl(Age_group, CPG)  %>% 
  adorn_percentages()

prop.trend.test(cpg_pts, all_pts) # Trend test 

example.melt <- reshape2::melt(freq_table, id.vars="Age_group")
example.melt$variable <- factor(example.melt$variable,levels=c("No","Yes"))
example.melt$Age_group<-factor(example.melt$Age_group,levels=c("TARGET","Under 30","30s","40s","50s","60s","70s","Over 80"),label=c("TARGET","19-29","30-39","40-49","50-59","60-69","70-79","Over 80"))


cpg_prevalent <- ggplot(example.melt, aes(x=Age_group, y=value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity",color='black',width=0.9) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = paste0(round(value*100,2),"%")), 
            position = position_stack(vjust = 0.5), size = 3.5)+scale_fill_manual(values=c("#d3d3d3","#1262b3"))+
  sm_hgrid(legends = TRUE)+guides(fill=guide_legend(title="pFGVs in CPGs"))+ylab("Percent")+
  xlab("Age group")+ggplot2::annotate("text",x=6.5,y=1.05,
                                      label=expression('Cochran-Armitage Trend Test'~italic(P)~'< 0.001'),
                                      size=3)
