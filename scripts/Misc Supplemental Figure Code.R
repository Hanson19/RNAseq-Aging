#Introduction: Misc Supplemental Figures

#Purpose: Show how supplemental Figures 2 and 7 were created

#Created: 8/28/2024, KMH
#Last Edited: 8/28/24, KMH

#Supplemental Figure 2 Correlation Dat and Survival####
library(tidyverse)
library(ggpubr)

# Master Directory/
#   deseq_results/
#     All Analyses/
#       Combo_Analysis_Clustered_Z_scores_k28.csv

z_scores <- read.csv("deseq_results/All Analyses/Combo_Analysis_Clustered_Z_scores_k28.csv", header = TRUE)
tp_unique <- unique(z_scores[4:5]) %>% arrange(day)
cor(tp_unique$day, tp_unique$survivorship)
# -0.9809378 #Negative because survivial starts at 100 and goes down to 0
cor(tp_unique$day, 1:15)
# [1] 0.9990591

day_surv_cor <- 
  tp_unique %>% ggplot(aes(x=day, y=survivorship))+geom_point()+
  geom_smooth(method = "lm", se=FALSE)+
  scale_y_reverse(breaks = c(1,0.8,0.6,0.4,0.2,0), limits=c(1.1,0),labels=c("100%", "80%", "60%", "40%", "20%", "0%"))+
  scale_x_continuous(n.breaks = 6)+
  stat_cor(method = "pearson")+
  ylab("Survivorship\n(Percent Alive)")+xlab("Day")+
  theme_bw()

ggsave("Plots/SupFig2_day_surv_Cor.png", day_surv_cor, width = 6, height = 5, units = "in")

#Supplemental Figure 7 RNA Concentraion####
library(tidyverse)
library(ggpubr)

#Will need Supplementary Table B
sup_table_B <- read.table("supp_table_B.txt", header = TRUE)

sup_table_B_subset <- sup_table_B %>% filter(Used_for_library == "Yes") 
fly_age <- sup_table_B_subset[[3]]
rna_conc <- sup_table_B_subset[[6]]

cor.test(fly_age, rna_conc)
# Pearson's product-moment correlation
# 
# data:  fly_age and rna_conc
# t = -4.4906, df = 46, p-value = 4.747e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7228171 -0.3177756
# sample estimates:
#        cor 
# -0.5520639

rna_age_df <- cbind(fly_age, rna_conc)
rna_age_df <- as.data.frame(rna_age_df)
head(rna_age_df)

cor_con_day_plot <- rna_age_df %>% ggplot(aes(x=fly_age, y=rna_conc))+
  geom_point(size=3)+
  geom_smooth(method = "lm", se=FALSE, color="red")+
  xlab("Fly Age (Days)")+
  ylab("RNA Concentration (ng/ul)")+
  stat_cor(method = "pearson", label.x = 40, label.y.npc = "top")+
  theme_bw()

ggsave("Plots/SupFig7_RNA_Conc_Decreasing_age.png",cor_con_day_plot,width = 6, height = 5, units = "in")
