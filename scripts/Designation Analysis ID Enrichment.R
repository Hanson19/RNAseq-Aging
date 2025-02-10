#Introduction####
#Title: Analysis Enrichment in Designated Trajectories
#
#Purpose: Do we see specific clusters or specific trajectories enriched
#for genes that were identified in certain analyses. Specifically, are the complex
#genes more often identified in sampling point analyses
#
#Setup: This code will show how to make supplementary figures 4 and 5. Code
#requires going through Identifying_and_Clustering.R

#Created: 6/24/2024, KMH
#Last Edited: 8/09/2024, KMH

#Expression Designation and Analysis ID#####
#Sup Figure 4####
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

#Code expects a certain directory structure and files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       Sig_Genes_AnalysisID.csv
#       gene_color_id_ref.csv

#Read in data tables
analysis_ids <- read.csv("deseq_results/All Analyses/Sig_Genes_AnalysisID.csv", header = TRUE)
gene_clus_info <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)

#Add the Official ID for each gene to the genes Analysis ID information
genes_clus_ids <- merge(analysis_ids[2:6], gene_clus_info[c(3,7)], by=c("gene_id"))
genes_clus_ids$Designation <- gsub("-.*","", genes_clus_ids$Official_ID)

#Note regarding Analysis ID:
#D stands for Day, the gene was identified in the chronological aging analysis
#S stands for Survival, the gene was identified in the physiological aging analysis
#P stands for Sampling Point, the gene was identified in the categorical analysis

#A gene can have any combination of these three letters. For example:
# DP means the gene was identified in day and sampling point, whereas
# D means the gene was only identified in day.

#Count how many Genes within a designation has a certain Analysis ID
desg_id_count <- genes_clus_ids %>% dplyr::count(Designation, Analysis_ID)
#Get the total number of genes for each Designation and each Analysis ID
desg_tg <- genes_clus_ids %>% dplyr::count(Designation)
colnames(desg_tg)[2] <- "Total_Genes_Designation"
id_tg <- genes_clus_ids %>% dplyr::count(Analysis_ID)
colnames(id_tg)[2] <- "Total_Genes_Analysis_ID"
desg_id_count <- merge(desg_id_count, desg_tg,by=c("Designation"))
desg_id_count <- merge(desg_id_count, id_tg,by=c("Analysis_ID"))
#Get the percent makeup of Designation/Analysis ID for Designation and Anlysis ID 
desg_id_count$Percent_Designation <- desg_id_count$n/desg_id_count$Total_Genes_Designation*100
desg_id_count$Percent_Analysis_ID <- desg_id_count$n/desg_id_count$Total_Genes_Analysis_ID*100

desg_id_count$Analysis_ID <- factor(desg_id_count$Analysis_ID, levels = c("D", "S", "P",
                                                                          "DP", "SP", "SD",
                                                                          "SDP"))

#Plotting the makeup for each Analysis ID
analysisid_desg<- 
  desg_id_count %>% ggplot(aes(x=Analysis_ID, y=Percent_Analysis_ID, fill=Designation))+
  geom_bar(stat = "identity")+
  scale_x_discrete(labels = c("Day", "Survival", "Sampling Point",
                              "Day\nSampling Point", "Survival\nSampling Point", "Day\nSurvival",
                              "Day\nSurvival\nSampling Point"))+
  xlab("Analyses Gene is Identified in")+
  scale_y_continuous(labels = c("0%", "25%", "50%", "75%","100%"))+
  scale_fill_manual(values = brewer.pal(11, "Paired")[c(4,2,6)], name="Trajectory Designation")+
  ylab("Percent")+
  ggtitle("A.")+
  theme_bw()+
  theme(text = element_text(size=7))

#Plotting the makeup for each Expression Designation
desg_analysisid<- 
  desg_id_count %>% ggplot(aes(x=Designation, y=Percent_Designation, fill=Analysis_ID))+
  geom_bar(stat = "identity",position = position_dodge(0.9))+
  geom_text(aes(label=n, y=Percent_Designation+1), position = position_dodge(0.9), size=2)+
  scale_fill_manual(values = brewer.pal(8,"Set2")[c(1:6,8)],labels=c("Day", "Survival", "Sampling Point",
                                                                     "Day+Sampling Point", "Survival+Sampling Point", "Day+Survival",
                                                                     "Day+Survival+\nSampling Point"),
                    name="Analysis")+
  scale_y_continuous(labels = c("0%", "20%", "40%", "60%","X"))+
  ylab("Percent")+
  xlab("Trajectory Designation")+
  ggtitle("B.")+
  theme_bw()+
  theme(text = element_text(size=7))

#As we can see in these two graphs the majority of the genes found in complex designations
#are genes that are only found in Sampling Point Analysis. Additionally we see that for
#genes identified in either just the sampling point analysis or samplingpoint with one
#of the other analyses the majority of the genes are found in complex clusters. 

#Combine and Save
desg_analysis_combine <- grid.arrange(analysisid_desg, desg_analysisid, ncol=1)
ggsave("Plots/SupFig4_Desg_Analysis_ID.png",desg_analysis_combine, width = 6.5, height = 6, units = "in")


#Cluster and Analysis ID####
#Sup Figure 5 Indv Clusters####
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

#Code expects a certain directory structure and files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       Sig_Genes_AnalysisID.csv
#       gene_color_id_ref.csv

#Read in data tables
analysis_ids <- read.csv("deseq_results/All Analyses/Sig_Genes_AnalysisID.csv", header = TRUE)
gene_clus_info <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)

#Add the Official ID for each gene to the genes Analysis ID information
genes_clus_ids <- merge(analysis_ids[2:6], gene_clus_info[c(3,7)], by=c("gene_id"))
genes_clus_ids$Designation <- gsub("-.*","", genes_clus_ids$Official_ID)

#Count how many genes within each cluster have specific Analysis ID
clus_id_count <- genes_clus_ids %>% dplyr::count(Official_ID, Analysis_ID,Designation)
#Get the Total number of genes for each cluster and each analysis ID
clus_tg <- genes_clus_ids %>% dplyr::count(Official_ID)
colnames(clus_tg)[2] <- "Total_Genes_Cluster"
id_tg <- genes_clus_ids %>% dplyr::count(Analysis_ID)
colnames(id_tg)[2] <- "Total_Genes_Analysis_ID"
clus_id_count <- merge(clus_id_count, clus_tg,by=c("Official_ID"))
clus_id_count <- merge(clus_id_count, id_tg,by=c("Analysis_ID"))
#Get the percent makeup of Cluster/Analysis ID for Cluster and Anlysis ID 
clus_id_count$Percent_Cluster <- clus_id_count$n/clus_id_count$Total_Genes_Cluster*100
clus_id_count$Percent_Analysis_ID <- clus_id_count$n/clus_id_count$Total_Genes_Analysis_ID*100

clus_id_count$Analysis_ID <- factor(clus_id_count$Analysis_ID, levels = c("D", "S", "P",
                                                                          "DP", "SP", "SD",
                                                                          "SDP"))


clus_id_count$Official_ID <- factor(clus_id_count$Official_ID, levels  = c("Complex-1","Complex-2","Complex-3","Complex-4",
                                                                       "Complex-5","Complex-6","Complex-7","Complex-8",
                                                                       "Complex-9", "Complex-10", "Complex-11","Complex-12",
                                                                       "Complex-13","Complex-14", "LinearUp-1","LinearUp-2",
                                                                       "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                                                       "LinearUp-7", "LinearUp-8", "LinearDown-1", "LinearDown-2",
                                                                       "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6"))


clus_id_count$color[clus_id_count$Designation == "Complex"] <- "white"
clus_id_count$color[clus_id_count$Designation == "LinearDown"] <- "gray65"
clus_id_count$color[clus_id_count$Designation == "LinearUp"] <- "gray90"

#Ploting the Analysis ID makeup for each cluster
cluster_analysis_id <- 
  clus_id_count %>% ggplot(aes(x=Analysis_ID, y=Percent_Cluster))+
  facet_wrap(~Official_ID,ncol = 7)+
  geom_bar(aes(fill=Analysis_ID),stat = "identity",position = position_dodge(0.9), show.legend = TRUE)+
  scale_fill_manual(values = brewer.pal(8,"Set2")[c(1:6,8)],labels=c("Day", "Survival", "Sampling Point",
                                                                     "Day+Sampling Point", "Survival+Sampling Point", "Day+Survival",
                                                                     "Day+Survival+\nSampling Point"),
                    name="Analysis")+
  scale_y_continuous(labels = c("0%", "25%", "50%", "75%","100%"))+
  ylab("Percent")+
  theme_bw()+
  xlab("Cluster")+
  coord_cartesian(clip="off", ylim=c(0, 100)) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=105, ymax=135),alpha=0.4,fill=clus_id_count$color)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill=NA),
        legend.text = element_text(size=7),
        legend.title = element_text(size=9),
        legend.position = "bottom"
        #text = element_text(size=10)
        )

#In this plot we can see the enrichment of sampling point only genes in complex
#clusters is not driven by one or two clusters, but is seen in all but 2 complex cluster.
#In addition we can see for all the linear cluster their top two analysis IDs for all 
#clusters are genes found in all three analyses or genes found only in Day and Survival. 

#Save plot
ggsave("Plots/SupFig5_Cluster_Analysis_ID.png",cluster_analysis_id, width = 6.5, height = 6, units = "in")
