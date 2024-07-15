# 1. Library ----------------------------------------------------------------------------------
package_list <- c("readr", 'tidyverse',"tibble", "openxlsx", "readxl", 
                  "writexl", "stringr", "riskRegression", "MASS", "devtools", 
                  "survAUC", "cowplot", "ggpubr", "survival", "survminer",
                  "rms", "prodlim", "RColorBrewer", "scales", "lava",'smplot2')
purrr::walk(package_list, library, character.only = TRUE)

output_folder <- '/path/output/SMC'
data_folder <-'/path/data/SMC/'
# 2. Data -------------------------------------------------------------------------------------

SMC_data <- readRDS(paste0(data_folder,'/SMC_data.Rds'))

# 3. Training ---------------------------------------------------------------------------------

y <- SMC_data

set.seed(1235)
iter <- 500
error_mat <- matrix(NaN, nrow=(iter*6), ncol=2)
colnames(error_mat) <- c("group", "cindex")
rownames(error_mat) <- c(1:(iter*6))

j <- 1

for(i in (1:iter)) {
  
  # permutation
  t <- 1:dim(y)[1]
  inx <- sample(t)
  
  x <- y[inx,]
  
  TR <- x[1:round(dim(x)[1]*0.60),]
  TE <- x[(round(dim(x)[1]*0.60)+1):dim(x)[1],]
  
  # Clinical
  try({fit <- rms::cph( Surv(fu_event_pr,event_pr=="Yes") ~ Age2+mycn+stage, data=TR, x=T, y=T,)
  val <- rms::validate( fit, dxy=T, B=200, )
  c.index1 <- 0.5 * ( val[1,1] + 1 ) 
  error_mat[j,1] <- 1
  error_mat[j,2] <- c.index1},silent=TRUE)
  j <- j+1
  
  # Germline variant burden
  try({fit <- rms::cph( Surv(fu_event_pr,event_pr=="Yes") ~ germ_burden_group, data=TR, x=T, y=T )
  val <- rms::validate( fit, dxy=T, B=200 )
  c.index2 <- 0.5 * ( val[1,1] + 1 ) 
  error_mat[j,1] <- 2
  error_mat[j,2] <- c.index2},silent=TRUE)
  j <- j+1
  
  # pFGV in CPGs  
  try({fit <- rms::cph( Surv(fu_event_pr,event_pr=="Yes") ~ CPG, data=TR, x=T, y=T )
  val <- rms::validate( fit, dxy=T, B=200 )
  c.index3 <- 0.5 * ( val[1,1] + 1 ) 
  error_mat[j,1] <- 3
  error_mat[j,2] <- c.index3},silent=TRUE)
  j <- j+1
  
  # CPG + Germline variant burden
  try({fit <- rms::cph( Surv(fu_event_pr,event_pr=="Yes") ~ germ_burden_group + CPG, data=TR, x=T, y=T )
  val <- rms::validate( fit, dxy=T, B=200 )
  c.index4 <- 0.5 * ( val[1,1] + 1 ) 
  error_mat[j,1] <- 4
  error_mat[j,2] <- c.index4},silent=TRUE)
  j <- j+1
  
  # Clinical + Germline variant burden + CPG 
  try({fit <- cph( Surv(fu_event_pr,event_pr=="Yes") ~Age2+mycn+stage + germ_burden_group +CPG, data=TR, x=T, y=T )
  val <- rms::validate( fit, dxy=T, B=200 )
  c.index5 <- 0.5 * ( val[1,1] + 1 )
  error_mat[j,1] <- 5
  error_mat[j,2] <- c.index5},silent=TRUE)
  j <- j+1
}


# 3. Testing ----------------------------------------------------------------------------------

set.seed(1235)
iter <- 500
error_mat_te <- matrix(NaN, nrow=(iter*6), ncol=2)
colnames(error_mat_te) <- c("group", "cindex")
rownames(error_mat_te) <- c(1:(iter*6))

j <- 1

for(i in (1:iter)) {
  
  # permutation
  t <- 1:dim(y)[1]
  inx <- sample(t)
  
  x <- y[inx,]
  
  TR <- x[1:round(dim(x)[1]*0.60),]
  TE <- x[(round(dim(x)[1]*0.60)+1):dim(x)[1],]
  
  # Clinical
  try({fit <- rms::cph( Surv(fu_event_pr,event_pr=="Yes") ~ Age2+mycn+stage, data=TE, x=T, y=T,)
  val <- rms::validate( fit, dxy=T, B=200, )
  c.index1 <- 0.5 * ( val[1,1] + 1 ) 
  error_mat_te[j,1] <- 1
  error_mat_te[j,2] <- c.index1},silent=TRUE)
  j <- j+1
  
  # Germline variant burden
  try({fit <- rms::cph( Surv(fu_event_pr,event_pr=="Yes") ~ germ_burden_group, data=TE, x=T, y=T )
  val <- rms::validate( fit, dxy=T, B=200 )
  c.index2 <- 0.5 * ( val[1,1] + 1 ) 
  error_mat_te[j,1] <- 2
  error_mat_te[j,2] <- c.index2},silent=TRUE)
  j <- j+1
  
  # pFGV in CPGs  
  try({fit <- rms::cph( Surv(fu_event_pr,event_pr=="Yes") ~ CPG, data=TE, x=T, y=T )
  val <- rms::validate( fit, dxy=T, B=200 )
  c.index3 <- 0.5 * ( val[1,1] + 1 ) 
  error_mat_te[j,1] <- 3
  error_mat_te[j,2] <- c.index3},silent=TRUE)
  j <- j+1
  
  # CPG + Germline variant burden
  try({fit <- rms::cph( Surv(fu_event_pr,event_pr=="Yes") ~ germ_burden_group + CPG, data=TE, x=T, y=T )
  val <- rms::validate( fit, dxy=T, B=200 )
  c.index4 <- 0.5 * ( val[1,1] + 1 ) 
  error_mat_te[j,1] <- 4
  error_mat_te[j,2] <- c.index4},silent=TRUE)
  j <- j+1
  
  # Clinical + Germline variant burden + CPG 
  try({fit <- cph( Surv(fu_event_pr,event_pr=="Yes") ~Age2+mycn+stage + germ_burden_group +CPG, data=TE, x=T, y=T )
  val <- rms::validate( fit, dxy=T, B=200 )
  c.index5 <- 0.5 * ( val[1,1] + 1 )
  error_mat_te[j,1] <- 5
  error_mat_te[j,2] <- c.index5},silent=TRUE)
  j <- j+1
}


# 4. Training Figure --------------------------------------------------------------------------
e_mat_tr["group"][e_mat_tr["group"] =="1"] <- "Clinical risk factors"
e_mat_tr$group[e_mat_tr$group == '2'] <- 'Germline risk factors (Burden)'
e_mat_tr$group[e_mat_tr$group == '3'] <- 'Germline risk factors (CPG)'
e_mat_tr$group[e_mat_tr$group == '4'] <- 'Germline risk factors (All)'
e_mat_tr$group[e_mat_tr$group == '5'] <- 'Clinical risk factors + Germline risk factors (All)'
e_mat_tr <-e_mat_tr %>% filter(!is.na(cindex))
e_mat_tr$group <- factor(e_mat_tr$group,levels=c("Germline risk factors (Burden)","Germline risk factors (CPG)",
                                                 "Germline risk factors (All)","Clinical risk factors",
                                                 "Clinical risk factors + Germline risk factors (All)"))

my_comparisons=list(c("Germline risk factors (All)","Clinical risk factors"),c("Clinical risk factors","Clinical risk factors + Germline risk factors (All)"),c("Clinical risk factors + Germline risk factors (All)","Germline risk factors (All)"))

my_comparisons=list(c("Germline risk factors (All)","Clinical risk factors"),
                    c("Clinical risk factors + Germline risk factors (All)","Germline risk factors (All)"),c("Clinical risk factors","Clinical risk factors + Germline risk factors (All)"))

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


c_index_tr <- e_mat_tr %>% 
  ggplot(aes(x=group, fill= group, y=cindex))+
  geom_point(position=position_jitterdodge(),
             color="black",
             shape=21,
             # fill="gray",
             stroke=0.2)+
  sm_boxplot(point.params = list(shape = NA, color = 'white',
                                 size = 1))+
  scale_fill_manual(values =sm_color("wine","skyblue","blue","green","viridian"))+
  ggpubr::stat_compare_means(method = "t.test",
                       
                             p.adjust.method = "bonferroni",
                             comparisons = my_comparisons,
                             label = "p.value",
  ) +

  theme_bw() +theme(legend.position = "NA",axis.title.x = element_text(),
                    axis.title = element_text(),axis.ticks.x =element_text(vjust=0.3),
                    axis.line = element_line(size=0.5, colour = "black"))+sm_hgrid()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size =10))+
  scale_x_discrete(breaks=unique(e_mat_tr$group), 
                   labels=c("Clinical risk factors","Germline risk factor\n(Burden)","Germline risk factor\n(CPG)",
                            "Germline risk factor\n(All)", 
                            "Clinical risk factors\n+\nGermline risk factors (All)"))+ylab("C-index")+xlab("")+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), expand = expansion(mult = c(0, .3)))


# 5. Testing ----------------------------------------------------------------------------------
e_mat_te["group"][e_mat_te["group"] =="1"] <- "Clinical risk factors"
e_mat_te$group[e_mat_te$group == '2'] <- 'Germline risk factors (Burden)'
e_mat_te$group[e_mat_te$group == '3'] <- 'Germline risk factors (CPG)'
e_mat_te$group[e_mat_te$group == '4'] <- 'Germline risk factors (All)'
e_mat_te$group[e_mat_te$group == '5'] <- 'Clinical risk factors + Germline risk factors (All)'
e_mat_te <-e_mat_te %>% filter(!is.na(cindex))
e_mat_te$group <- factor(e_mat_te$group,levels=c("Germline risk factors (Burden)","Germline risk factors (CPG)",
                                                 "Germline risk factors (All)","Clinical risk factors",
                                                 "Clinical risk factors + Germline risk factors (All)"))




c_index_te <- e_mat_te %>% 
  ggplot(aes(x=group, fill= group, y=cindex))+
  geom_point(position=position_jitterdodge(),
             color="black",
             shape=21,
             # fill="gray",
             stroke=0.2)+
  sm_boxplot(point.params = list(shape = NA, color = 'white',
                                 size = 1))+
  scale_fill_manual(values =sm_color("wine","skyblue","blue","green","viridian"))+
  ggpubr::stat_compare_means(method = "t.test",
                             p.adjust.method = "bonferroni",
                             comparisons = my_comparisons,
                             label = "p.value",
  ) +

  theme_bw() +theme(legend.position = "NA",axis.title.x = element_text(),
                    axis.title = element_text(),axis.ticks.x =element_text(vjust=0.3),
                    axis.line = element_line(size=0.5, colour = "black"))+sm_hgrid()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size =10))+
  scale_x_discrete(breaks=unique(e_mat_tr$group), 
                   labels=c("Clinical risk factors","Germline risk factor\n(Burden)","Germline risk factor\n(CPG)",
                            "Germline risk factor\n(All)", 
                            "Clinical risk factors\n+\nGermline risk factors (All)"))+ylab("C-index")+xlab("")+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), expand = expansion(mult = c(0, .3)))

