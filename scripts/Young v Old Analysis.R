#Introduction####
#Title: Simulated Old vs Young Analysis

#Purpose: To understand the impact of taking samples at multiple time points
#vs the traditional young vs old approach we are going to simulate a young
#vs old analysis with our samples. Young will be composed of samples taken
#at days 3 and 6 and old with samples from Day59.We look to see what genes are 
#shared between Old v Young and our original multiple time point analysis, if
#we see the genes being shared enriched for onr of our designated trajectories,
#and if the direction called in Old v Young analysis matches with out data. 
#In addition script shows how to make Figure 4 and Supplemental Figure 8

#Setup:
#Each section begins with brief description of the goal, what R packages are 
#required and the directory structure including files the code expects. Assumes 
#you have made gene_color_id_ref.csv (See Identifying_and_Clustering.R) and
#gse_salmon_tximeta_ovy.RDS

#Created:6/14/2024 KMH
#Last Edited: 8/28/2024 KMH

#Identifying Genes with DESeq####
library(BiocManager)
library(DESeq2)
library(tidyverse)
library(Rmisc)

#Code expects a certain directory structure with this specific files:
# Master Directory/
#   Aging Analysis Results/
#     Old v Young Sim/
#       gse_salmon_tximeta_20240521_ovy.RDS

#Read in gse
gse <- readRDS("Aging Analysis Results/Old v Young Sim/gse_salmon_tximeta_ovy.RDS")

str(metadata(gse))
colData(gse)
rowData(gse)

#Use DESeq to identify genes whose expression changes with age which is either
#young or old
dds <- DESeqDataSet(gse, design = ~ age)
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds)
sig_test <- res
res$symbol <- mcols(dds)$symbol
res$gene_id <- mcols(dds)$gene_id
sig_test$symbol <- mcols(dds)$symbol
sig_test$gene_id <- mcols(dds)$gene_id
#Filter so significant genes who padj is less than 0.05
sig_test <- sig_test[sig_test$padj < 0.05 &
                       !is.na(sig_test$pvalue) &
                       !is.na(sig_test$padj),]
sig_test <- sig_test[order(sig_test$padj),]
rownames(sig_test) <- 1:nrow(sig_test)

#Save all genes deseq results, and significant genes specific results
write.csv(as.data.frame(sig_test),"Aging Analysis Results/Old v Young Sim/OvY_PadjSig_deseq.csv")
write.csv(as.data.frame(res), "Aging Analysis Results/Old v Young Sim/OvY_All_deseq.csv")

#Save normalized counts
read_counts <- counts(dds, normalized=TRUE)
write.csv(read_counts,"Aging Analysis Results/Old v Young Sim/OvY_NormCounts_deseq.csv")

dim(sig_test)
# [1] 4533    8

#Shared Significant genes####
#How many genes are shared between our multi time point analysis and old v young
library(tidyverse)

#Expected Directory structure and files
# Master Directory/
#   Aging Analysis Results/
#     Old v Young Sim/
#       OvY_PadjSig_deseq.csv
#   deseq_results/
#     All Analyses/
#       gene_color_id_ref.csv

YvO_DE_Genes <- read.csv("Aging Analysis Results/Old v Young Sim/OvY_PadjSig_deseq.csv", header = TRUE)
#Total: 4533 genes
Aging_DE_Genes_IDs <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)
#Total: 6142

#3 situations we could have
#1. Gene is found both in time analysis and in YvO
#2. Gene is unique to time analysis
#3. Gene is unique to YvO analysis

#Genes shared by both analyses:
shared_genes <- Aging_DE_Genes_IDs %>% filter(gene_id %in% YvO_DE_Genes$gene_id)
#4347 genes shared

#Genes unique to multi-time point analysis
time_genes <- Aging_DE_Genes_IDs %>% filter(!gene_id %in% YvO_DE_Genes$gene_id)
#1795 genes (~30% of our identified genes)

#Genes unique to YvO analysis
YvO_genes <- YvO_DE_Genes %>% filter(!gene_id %in% Aging_DE_Genes_IDs$gene_id)
#186 genes (~4% of YvO identified genes)


#Designated Trajectory Enrichment####
#So yes there are differences in genes being identified. We lost about 
#30% of our genes, so question now is are those lost genes found more
#concentrated in one of our Designated trajectories.
#Here we will compare Complex trajectories vs our linear trajectories
library(tidyverse)

#Expected Directory structure and files
# Master Directory/
#   Aging Analysis Results/
#     Old v Young Sim/
#       OvY_PadjSig_deseq.csv
#   deseq_results/
#     All Analyses/
#       gene_color_id_ref.csv

#Read in OvY signficiant genes and multi timepoint significant genes
YvO_DE_Genes <- read.csv("Aging Analysis Results/Old v Young Sim/OvY_PadjSig_deseq.csv", header = TRUE)
#Total: 4533 genes
Aging_DE_Genes_IDs <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)
#Total: 6142

#Filter out the unique time point genes
time_genes <- Aging_DE_Genes_IDs %>% filter(!gene_id %in% YvO_DE_Genes$gene_id)
time_genes$Official_ID <- factor(time_genes$Official_ID, levels = c("Complex-1","Complex-2","Complex-3","Complex-4","Complex-5","Complex-6",
                                                                    "LinearUp-1",
                                                                    "Complex-7",
                                                                    "LinearUp-2","LinearUp-3","LinearUp-4", "LinearUp-5",
                                                                    "Complex-8",
                                                                    "LinearUp-6", "LinearUp-7","LinearUp-8",
                                                                    "Complex-9", "Complex-10","Complex-11",
                                                                    "LinearDown-1", "LinearDown-2", "LinearDown-3", "LinearDown-4","LinearDown-5",
                                                                    "Complex-12",
                                                                    "LinearDown-6",
                                                                    "Complex-13","Complex-14"))
#Get total gene count for each cluster
TG <- Aging_DE_Genes_IDs %>% dplyr::count(Official_ID)
#Get the total unique gene for each cluster
UG <- time_genes %>% dplyr::count(Official_ID)
#Merge together
cluster_unique_percent <- merge(TG, UG, by=c("Official_ID"))
colnames(cluster_unique_percent) <- c("Official_ID", "Total_Genes", "Unique_Genes")
#Get a cluster's percent makeup of unique genes
cluster_unique_percent$Percent <- cluster_unique_percent$Unique_Genes/cluster_unique_percent$Total_Genes
#Get a cluster's designation of being Complex or Linear by removing the number from the name
cluster_unique_percent$Designation <- gsub("-.*","", cluster_unique_percent$Official_ID)
#Make LinearUp and LinearDown just linear by identifying cluster's whose designation beings with
#"Linear"
for (i in 1:nrow(cluster_unique_percent)) {
  if(str_starts(cluster_unique_percent[i,5], "Linear")){
    cluster_unique_percent[i,5] <- "Linear"
  }
}
#T.test comparing percent makeup of Linear Clusters and Complex clusters
t.test(cluster_unique_percent %>% filter(Designation == "Linear") %>% pull(Percent),
       cluster_unique_percent %>% filter(Designation == "Complex") %>% pull(Percent))
# p-value = 0.001787 There is a significant difference in makeup of unique genes for 
# clusters that are complex or linear


#Plot each Cluster Makeup
cluster_unique_percent %>% ggplot(aes(x=Designation, y= Percent))+
  geom_boxplot()+
  ylab("Proportion Unique")+
  geom_text(aes(label=Official_ID),position = position_jitter(.25))+
  theme_classic()
#Can see that for majority of out complex clusters .5 of its makeup is made
#of unique genes. Complex 12, 9, 6 and 7 are clusters that when looking at
#rep curves are all curves where there is a substantial difference in where 
#expression started vs ended. Complex 9, 6, 7 are also clusters whose
#linear regression with aging p-value was below 0.05, but above 0.002

#Fig 4 Designated Trajectory Direction####
#Old v Young anlaysis can only call if a gene's expression increases
#or decreases. Our question is how do these call stack up within our
#designated trajectories. Do we find the genes in our linearup trajectories, genes
#whose expression was called as increasing. 

library(tidyverse)

#Are the directions called in YvO analysis agree with our trajectories
#Expected Directory structure and files
# Master Directory/
#   Aging Analysis Results/
#     Old v Young Sim/
#       OvY_PadjSig_deseq.csv
#       OvY_NormCounts_deseq.csv
#   deseq_results/
#     All Analyses/
#       gene_color_id_ref.csv

#Read in OvY Padj Sig genes, normalized counts, and multi time point siginfiicant gene's ids
YvO_DE_Genes <- read.csv("Aging Analysis Results/Old v Young Sim/OvY_PadjSig_deseq.csv", header = TRUE)
YvO_Counts <- read.csv("Aging Analysis Results/Old v Young Sim/OvY_NormCounts_deseq.csv", header = TRUE)
Time_Sig_Genes <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)

#Get the normalized counts for only the OvY significant genes
Sig_YVO_Counts <- YvO_Counts %>% filter(X %in% YvO_DE_Genes$gene_id)
colnames(Sig_YVO_Counts)[1] <- "gene_id"
Sig_YVO_Counts$Young_Count <- NA
Sig_YVO_Counts$Old_Count <- NA
Sig_YVO_Counts$Direction <- NA
#For each gene (Each gene is a row)
for (r in 1:nrow(Sig_YVO_Counts)) {
  #Get the mean count for our young counts
  Sig_YVO_Counts[r,14] <- mean(Sig_YVO_Counts[r,2],Sig_YVO_Counts[r,3],
                               Sig_YVO_Counts[r,4],Sig_YVO_Counts[r,5],
                               Sig_YVO_Counts[r,6],Sig_YVO_Counts[r,7])
  #Get the mean count for our old counts
  Sig_YVO_Counts[r,15] <- mean(Sig_YVO_Counts[r,8],Sig_YVO_Counts[r,9],
                               Sig_YVO_Counts[r,10],Sig_YVO_Counts[r,11],
                               Sig_YVO_Counts[r,12],Sig_YVO_Counts[r,13])
  #If young mean is greater than old mean that call the direction down
  if(Sig_YVO_Counts[r,14] > Sig_YVO_Counts[r,15]){
    Sig_YVO_Counts[r,16] <- "Down"
  }else{
  #Else call the direction Up 
    Sig_YVO_Counts[r,16] <- "Up"
  }
}

write.csv(Sig_YVO_Counts, "Aging Analysis Results/Old v Young Sim/Sig_Gene_Direction_Change_59v36.csv")

#Merge counts with multi time point significant gene info and keep the
#young and old counts, the called direction, and cluster info and analysis ID from
#multi time point analysis
Sig_YVO_Counts <- merge(Sig_YVO_Counts, Time_Sig_Genes[2:8], by="gene_id")
Sig_YVO_Counts<- Sig_YVO_Counts[,c(1,14:19,21)]

#Get the cluster's Designation by removing the number from official ID
Sig_YVO_Counts$Designation <- gsub("-.*","", Sig_YVO_Counts$Official_ID)

#count the number of genes within each classified trajectory that was called up or down
desg_direction_sum <- Sig_YVO_Counts %>% dplyr::count(Designation, Direction)
#   Designation Direction    n
# 1     Complex      Down  222
# 2     Complex        Up  288
# 3  LinearDown      Down 1987
# 4  LinearDown        Up   24
# 5    LinearUp      Down   33
# 6    LinearUp        Up 1793

#Get number of genes not identified in the old v young analysis
#Pull out genes not identified in ovy analysis
Not_shared <- Time_Sig_Genes %>% filter(!gene_id %in% Sig_YVO_Counts$gene_id)
#Get cluster's designation and then count how many genes for each designation.
Not_shared$Designation <- gsub("-.*","", Not_shared$Official_ID)
Not_shared_sum <- Not_shared %>% dplyr::count(Designation)
#Give it a direction of "Not identified
Not_shared_sum$Direction <- "Not Identified"

#Combine Identified and not identified gene counts
desg_direction_sum <- rbind(desg_direction_sum, Not_shared_sum)

#Get total number of genes within each designation
Time_Sig_Genes$Designation <- gsub("-.*","", Time_Sig_Genes$Official_ID)
TG_Count <- Time_Sig_Genes %>% dplyr::count(Designation)
colnames(TG_Count)[2] <- "Total_Genes"

#Add Total number of genes to desg_direction_sum
desg_direction_sum <- merge(desg_direction_sum, TG_Count, by="Designation")
#Get the percent makeup of each direction for that designation's gene pool
desg_direction_sum$Percent <- paste(round(desg_direction_sum$n/desg_direction_sum$Total_Genes*100,2), "%", sep = "")

#Plot the makeup of each trajectory designation 
desg_direction_sum$Label_place <- c(111,366,760,
                                    1999,994,2304,
                                    17,930,2181)
desg_direction_sum$Direction <- factor(desg_direction_sum$Direction, levels = c("Not Identified",
                                                                                "Up",
                                                                                "Down"))

trajectory_direction <- 
  desg_direction_sum %>% ggplot(aes(x=Designation, y=n, fill=Direction))+
  geom_bar(stat = "identity")+
  xlab("Cluster Expression Trajectory")+
  ylab("Number of Genes")+
  scale_fill_manual(values = c("gray93", "skyblue", "tomato"))+
  geom_label(aes(label=Percent, y=Label_place), show.legend = FALSE, 
             label.size = 0, fill="white", alpha=0, size=3)+
  guides(fill=guide_legend(title = c("Old v Young\nDirection Call")))+
  ggtitle("Old v Young Expression Change Designation")+
  theme_light()+
  theme(text = element_text(size = 10))

#When looking at graph we can see for our linear designations the called direction
#in old v young analysis matches the direction of linear designation for all but 
#about 1% of that designation's genes (its between 1% and 2% when looking at just
#shared genes). For our Complex genes 50% of the genes in that pool are not
#identified young v old analysis and for the genes that are identified they are
#heavily split between being called Up or Down. 

ggsave("Plots/Fig4_Trajectory_ovy_direction.png", trajectory_direction, height = 6, width = 6, units = "in")


#Sup Fig 8 YvO Breakdown Indv Clusters####
library(tidyverse)
#Plot for each cluster the percent composition of genes called up and down in 
#young v old analysis and those not identified. 

#Expected Directory structure and files
# Master Directory/
#   Aging Analysis Results/
#     Old v Young Sim/
#       OvY_PadjSig_deseq.csv
#       OvY_NormCounts_deseq.csv
#   deseq_results/
#     All Analyses/
#       gene_color_id_ref.csv

#Read in OvY Padj Sig genes, normalized counts, and multi time point siginfiicant gene's ids
YvO_DE_Genes <- read.csv("Aging Analysis Results/Old v Young Sim/OvY_PadjSig_deseq.csv", header = TRUE)
YvO_Counts <- read.csv("Aging Analysis Results/Old v Young Sim/OvY_NormCounts_deseq.csv", header = TRUE)
Time_Sig_Genes <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)

#Get the normalized counts for only the OvY significant genes
Sig_YVO_Counts <- YvO_Counts %>% filter(X %in% YvO_DE_Genes$gene_id)
colnames(Sig_YVO_Counts)[1] <- "gene_id"
Sig_YVO_Counts$Young_Count <- NA
Sig_YVO_Counts$Old_Count <- NA
Sig_YVO_Counts$Direction <- NA
#For eahc gene (Each gene is a row)
for (r in 1:nrow(Sig_YVO_Counts)) {
  #Get the mean count for our young counts
  Sig_YVO_Counts[r,14] <- mean(Sig_YVO_Counts[r,2],Sig_YVO_Counts[r,3],
                               Sig_YVO_Counts[r,4],Sig_YVO_Counts[r,5],
                               Sig_YVO_Counts[r,6],Sig_YVO_Counts[r,7])
  #Get the mean count for our old counts
  Sig_YVO_Counts[r,15] <- mean(Sig_YVO_Counts[r,8],Sig_YVO_Counts[r,9],
                               Sig_YVO_Counts[r,10],Sig_YVO_Counts[r,11],
                               Sig_YVO_Counts[r,12],Sig_YVO_Counts[r,13])
  #If young mean is greater than old mean that call the direction down
  if(Sig_YVO_Counts[r,14] > Sig_YVO_Counts[r,15]){
    Sig_YVO_Counts[r,16] <- "Down"
  }else{
    #Else call the direction Up 
    Sig_YVO_Counts[r,16] <- "Up"
  }
}

#Merge counts with multi time point significant gene info and keep the
#young and old counts, the called direction, and cluster info and analysis ID from
#multi time point analysis
Sig_YVO_Counts <- merge(Sig_YVO_Counts, Time_Sig_Genes[2:8], by="gene_id")
Sig_YVO_Counts<- Sig_YVO_Counts[,c(1,14:19,21)]

#For each cluster count how many genes were called up and down and how many
#were not identified in young vs old analysis
clust_count <- Sig_YVO_Counts %>% dplyr::count(Official_ID, Direction)
clust_not <- Time_Sig_Genes %>% filter(!gene_id %in% Sig_YVO_Counts$gene_id)
clust_not_count <- clust_not %>% dplyr::count(Official_ID)
clust_not_count$Direction <- "Not Identified"
clust_count <- rbind(clust_count,clust_not_count)

#Get the total number flies per cluster
TG_Clust <- Time_Sig_Genes %>% dplyr::count(Official_ID)
clust_count <- merge(clust_count, TG_Clust, by=c("Official_ID"))
#Divide number of genes for each cluster's direction by a cluster's total number of
#genes then times by 100 to get the percent value
clust_count$Percent <- paste(round(clust_count$n.x/clust_count$n.y*100, 2), "%", sep = "")
clust_count$Percent_Num <- round(clust_count$n.x/clust_count$n.y*100,2)

#Plot Factorize IDs and direction so plotted in desired formation
clust_count$Direction <- factor(clust_count$Direction, levels = c("Not Identified","Up","Down"))
clust_count$Official_ID <- factor(clust_count$Official_ID, levels  = c("Complex-1","Complex-2","Complex-3","Complex-4",
                                                                                "Complex-5","Complex-6","Complex-7","Complex-8",
                                                                                "Complex-9", "Complex-10", "Complex-11","Complex-12",
                                                                                "Complex-13","Complex-14", "LinearUp-1","LinearUp-2",
                                                                                "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                                                                "LinearUp-7", "LinearUp-8", "LinearDown-1", "LinearDown-2",
                                                                                "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6"))
#Calculate label placement for plots
clust_count$label_place <- NA
for (c in unique(clust_count$Official_ID)) {
  d <- clust_count %>% filter(Official_ID == c & Direction == "Down") %>% pull(n.x)
  u <- clust_count %>% filter(Official_ID == c & Direction == "Up") %>% pull(n.x)
  ni <- clust_count %>% filter(Official_ID == c & Direction == "Not Identified") %>% pull(n.x)
  if(length(u) == 0){
    u <- 0
  }
  if(length(d) == 0){
    d <- 0
  }
  clust_count$label_place[clust_count$Official_ID == c & clust_count$Direction == "Down"] <- d/2
  clust_count$label_place[clust_count$Official_ID == c & clust_count$Direction == "Up"] <- d+u/2
  clust_count$label_place[clust_count$Official_ID == c & clust_count$Direction == "Not Identified"] <- d+u+ni/2
  # clust_count$label_place[clust_count$Official_ID == c & clust_count$Direction == "Down"] <- d/2
  # clust_count$label_place[clust_count$Official_ID == c & clust_count$Direction == "Up"] <- u/2
  # clust_count$label_place[clust_count$Official_ID == c & clust_count$Direction == "Not Identified"] <- ni/2
}

#Put placement in terms of percentage
clust_count$label_place_perc <- clust_count$label_place/clust_count$n.y*100

#Want to color Official Ids based on their designated trajectory. 
#Create a dataframe with that information
facet_df <- as.data.frame(clust_count$Official_ID)
facet_df$Designation <- gsub("-.*","", clust_count$Official_ID)
colnames(facet_df)[1] <- "Official_ID"

#Slight adjustment to LinearDown-6's label for Up genes, so labels
#are not on top of each other
clust_count[51,8] <- 82 #OG90.454545

#Plot barplot for each cluster showing their percent composition of OldvYoung calls
#have each section labeled with the number of genes in that group
cluster_breakdown_plot <- 
  clust_count %>% 
  ggplot()+
  #Color the facet wrap labels' background based on the cluster's trajectory
  geom_rect(data = facet_df, aes(xmin=-Inf, xmax=Inf, ymin=105, ymax=135,fill=Designation),alpha=0.4)+
  #Facet wrap based on official ID
  facet_wrap(~Official_ID, ncol = 7, scales = "free_x")+
  #Plot each cluster's percent composition 
  geom_bar(aes(x=Official_ID, y=Percent_Num, fill=Direction),stat = "identity", show.legend = TRUE)+
  xlab("Cluster")+
  ylab("Percent")+
  #labeling each section with the number of genes for each call wihtin cluster
  geom_label(aes(label=n.x, y=label_place_perc,x=Official_ID),
             label.size = 0, fill="white", alpha=0, size=3)+
  #Manually choosing colors and formatting legend so only the OldvYoung Call is
  #included in legend
  scale_fill_manual(values = c("white", #Complex
                                 "tomato", #Down
                                 "gray65", #Linear Down
                                 "gray90",  #Linear Up
                                 "gray87",  #Not identified
                                 "skyblue"), #Up
                      labels=c("Not Identified", "", "Up","","Down",""))+ #Up
  guides(fill=guide_legend(title = c("Young v. Old Analysis Call"),
                           title.vjust = .8,
                           override.aes = list(fill=c("gray87", "white", "skyblue", "white", "tomato", "white"))))+
  #Set limits of plot so can color facet wrap background
  coord_cartesian(clip="off", ylim=c(0, 100)) +
  theme_bw()+
  theme(text = element_text(size = 10), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = NA),
        legend.position = "bottom")

#Save plot
ggsave("Plots/SupFig8_IndvCluster_yvo.png",cluster_breakdown_plot, width = 6, height = 6, units = "in")
