#Introduction####
#Title: Identifying and Clustering Genes using Expression Trajectories

#Purpose: 
#This script will go through how genes were identified, how genes were 
#clustered together using expression trajectories, how number of clusters was determined
#and how a cluster's trajectory was designated. 
#In addition script shows how supplementary Figures 2-3, Figures 2-4 and 
#Supplementary Table D and E was made.

#Setup:
#Each section begins with what R packages are required and the directory
#structure including files the code expects. All files you will need will be
#made in the following sections except for the two files needed in "Gene Identification
#and Z Scores"

#Created: 5/30/2024, KMH
#Last Edited: 8/27/2024, KMH

#
#Gene Identification and Z Scores####
#Packages:
library(BiocManager)
library(DESeq2)
library(tidyverse)
library(Rmisc)
library(factoextra)
library(rms)

#Code expects a certain directory structure:
# Master Directory/
#   deseq_results/
#       day/
#       sampling_point/
#       survivorship/
#   Survivorship/
#       output.KaplanMeier.survivorship.txt #Code to create Aging_Survivorship_Code.R
#   Quant Files/
#       gse_salmon_tximeta.RDS

#Read in summarized Gene counts created using tximeta. 
gse <- readRDS("Quant Files/gse_salmon_tximeta.RDS")
str(metadata(gse))
colData(gse)
rowData(gse)

#We are identifying genes whose expression changes with aging with three different
#analyses:
#1. Survivorship: Identifying genes whose expression changes when against measured in terms of survivorship. Continuous Analysis
#2. Day: Identifying genes whose expression changes when measured against chronological time. Continuous Analysis
#3. Sampling Point:Identifying genes whose expression is variable among sampling point. Categorical Analysis

ID_and_Zscore <- function(GSE, ANALYSIS){
  #Setup DESeq Data Set with the design formula that is dependent on which Analysis is running
  if(ANALYSIS == "survivorship"){
    dds <- DESeqDataSet(GSE, design = ~ survivorship)
  }else if (ANALYSIS == "day"){
    dds <- DESeqDataSet(GSE, design = ~ day)
  }else if (ANALYSIS == "sampling_point"){
    gse$day <- factor(gse$day)
    dds <- DESeqDataSet(GSE, design = ~ day)
    dds$day <- factor(dds$day, levels = c("3","6","10","14","17",
                                          "23","27","31","36", 
                                          "38","42","48","50","55","59"))
  }
  dds <- DESeq(dds, test = "LRT", reduced = ~1)
  #pull out the results of DESeq Analysis
  res <- results(dds)
  sig_test <- res
  sig_test$symbol <- mcols(dds)$symbol
  sig_test$gene_id <- mcols(dds)$gene_id
  #pull out genes whose expression significantly changes using a threshold of padj < 0.05
  sig_test <- sig_test[sig_test$padj < 0.05 &
                         !is.na(sig_test$pvalue) &
                         !is.na(sig_test$padj),]
  sig_test <- sig_test[order(sig_test$padj),]
  rownames(sig_test) <- 1:nrow(sig_test)
  #Save both significant genes df, and reults df
  sig_test_filename <- paste("deseq_results/",ANALYSIS,"/",ANALYSIS,"_PadjSig_deseq.csv", sep = "")
  results_filename <- paste("deseq_results/",ANALYSIS,"/",ANALYSIS,"_All_deseq.csv", sep = "")
  write.csv(as.data.frame(sig_test), file = sig_test_filename)
  write.csv(as.data.frame(res), file = results_filename)
  
  #Now that genes have been identified we want to normalize gene expression using z scores 
  #so each gene's expression is on the same level
  #Pull out gene counts 
  read_counts <- counts(dds, normalized=TRUE)
  #The counts should be the same across the three analyses
  counts_filename <- paste("deseq_results/",ANALYSIS,"/",ANALYSIS,"_NormCounts_deseq.csv", sep = "")
  write.csv(read_counts, counts_filename)
  #Calculate Z_scores, only going to do this for our padj significant genes
  padj_read_counts <- as.data.frame(read_counts[sig_test$gene_id,])
  padj_read_counts <- rownames_to_column(padj_read_counts,"gene_id")
  read_counts_long <- padj_read_counts %>% pivot_longer(c(2:49), names_to = "Samples", 
                                                        values_to = "Norm_Counts")
  #Adding day and survivorship to table
  read_counts_long$day <- str_sub(read_counts_long$Samples, 1, -2)
  day_to_survivorship <- read.table("Survivorship/output.KaplanMeier.survivorship.txt", header = TRUE)
  included_day <- day_to_survivorship[day_to_survivorship$sampling.days == 1,1:2]
  colnames(included_day) <- c("day", "survivorship")
  read_counts_long <- merge(read_counts_long, included_day, by=c("day"))
  #Average each gene's read counts for each day
  gene_day_mean <- summarySE(data = read_counts_long, measurevar = c("Norm_Counts"), groupvars = c("gene_id","day","survivorship"))
  #Average each gene's read counts across every day
  gene_mean <- summarySE(data = gene_day_mean, measurevar = "Norm_Counts", groupvars = c("gene_id"))
  gene_Z_info <- gene_mean[c(1,3,4)]
  colnames(gene_Z_info) <- c("gene_id", "Mean_Counts_AllTP", "SD_AllTP")
  gene_Z_scores <- merge(gene_day_mean, gene_Z_info, by=c("gene_id"))
  #Subtract avg express across all days from the avg expression on one day and then divide by the standard diviation across all timepoints
  gene_Z_scores$Z_score <- (gene_Z_scores$Norm_Counts-gene_Z_scores$Mean_Counts_AllTP)/gene_Z_scores$SD_AllTP
  Z_scores_filename <- paste("deseq_results/",ANALYSIS,"/",ANALYSIS,"_Padj_Z_scores.csv", sep = "")
  write.csv(gene_Z_scores, file = Z_scores_filename)
}

A1 <- c("day")
A2 <- c("survivorship")
A3 <- c("sampling_point")
ID_and_Zscore(gse, A1)
ID_and_Zscore(gse, A2)
ID_and_Zscore(gse, A3)

#Should have 4 files in each of the analysis folder in deseq_results
# analysis_All_deseq
# analysis_NormCounts_deseq
# analysis_Padj_Z_score
# analysis_PadjSig_deseq

#Combined Analysis IDs####
#Create a document that has all the identified genes across all three analyses.
#As well as information regarding which analysis they were identified in. 
#Packages:
library(tidyverse)
library(ComplexUpset)

#Code expects a certain directory structure:
# Master Directory/
#   deseq_results/
#       day/
#         day_PadjSig_deseq.csv
#       sampling_point/
#         sampling_point_PadjSig_deseq.csv
#       survivorship/
#         survivorship_PadjSig_deseq.csv
#

survivorship <- read.csv("deseq_results/survivorship/survivorship_PadjSig_deseq.csv", header = TRUE)
sampling_point <- read.csv("deseq_results/sampling_point/sampling_point_PadjSig_deseq.csv", header = TRUE)
day <- read.csv("deseq_results/day/day_PadjSig_deseq.csv", header = TRUE)

#Create an ID system for genes to designate which time frames they were identified in
#Need to create an ID system
#Pull out gene ids for genes identified in the day and sampling point analysis
day_genes <- day[[9]]
sampling_points_genes <- sampling_point[[9]]

#first make list with unique genes
#Base list are genes identified in the survivorship analysis
all_sig_genes <- survivorship$gene_id
#add genes from the day analysis not found in the survivorship analysis
all_sig_genes <- append(all_sig_genes, day_genes[!day_genes %in% all_sig_genes])
#add genes from sampling point analysis that is unique to it. 
all_sig_genes <- append(all_sig_genes, sampling_point$gene_id[!sampling_point$gene_id %in% all_sig_genes])
all_sig_genes <- as.data.frame(all_sig_genes)
colnames(all_sig_genes) <- c("gene_id")

#Make an Analysis ID
#Mark if a gene was identified in a certain analysis 
all_sig_genes$Survivorship <- NA
all_sig_genes$Day <- NA
all_sig_genes$Sampling_Point <- NA

#If a gene was found in an anlysis it gets a letter in that analysis's column
all_sig_genes$Survivorship[all_sig_genes$gene_id %in% survivorship$gene_id] <- c("S") #S for survivorship
all_sig_genes$Day[all_sig_genes$gene_id %in% day$gene_id] <- c("D") #D for day
all_sig_genes$Sampling_Point[all_sig_genes$gene_id %in% sampling_point$gene_id] <- c("P") #P for sampling point

#If there is still NA in a column, replace the value with ""
all_sig_genes$Survivorship[is.na(all_sig_genes$Survivorship)] <- c("")
all_sig_genes$Day[is.na(all_sig_genes$Day)] <- c("")
all_sig_genes$Sampling_Point[is.na(all_sig_genes$Sampling_Point)] <- c("")

#Create an anlsysis ID by pasting together the 3 letters
all_sig_genes$Analysis_ID <- paste(all_sig_genes$Survivorship,
                                   all_sig_genes$Day, 
                                   all_sig_genes$Sampling_Point, sep = "")

dir.create("deseq_results/All Analyses")
#Save file
write.csv(all_sig_genes,"deseq_results/All Analyses/Sig_Genes_AnalysisID.csv")

##Figure 2 UpsetR of Analysis Identification####
#Make Upset figure to show the overlap of genes identified in different analyses
library(ComplexUpset)

Sig_Gene_Analysis_ID <- read.csv("deseq_results/All Analyses/Sig_Genes_AnalysisID.csv", header = TRUE,
                                 colClasses = "character")
Sig_Gene_Analysis_ID<- Sig_Gene_Analysis_ID[-c(1,6)]

#Create a table that for each of the analysis columns genes have 1s if they were
#identified in that anlysis and 0s if they were not

Sig_Gene_Analysis_ID$Sampling_Point[Sig_Gene_Analysis_ID$Sampling_Point  != "P"] <- 0
Sig_Gene_Analysis_ID$Sampling_Point[Sig_Gene_Analysis_ID$Sampling_Point == "P"] <- 1
Sig_Gene_Analysis_ID$Day[Sig_Gene_Analysis_ID$Day != "D"] <- 0
Sig_Gene_Analysis_ID$Day[Sig_Gene_Analysis_ID$Day == "D"] <- 1
Sig_Gene_Analysis_ID$Survivorship[Sig_Gene_Analysis_ID$Survivorship != "S"] <- 0
Sig_Gene_Analysis_ID$Survivorship[Sig_Gene_Analysis_ID$Survivorship == "S"] <- 1
Sig_Gene_Analysis_ID$Survivorship <- as.numeric(Sig_Gene_Analysis_ID$Survivorship)
Sig_Gene_Analysis_ID$Day <- as.numeric(Sig_Gene_Analysis_ID$Day)
Sig_Gene_Analysis_ID$Sampling_Point <- as.numeric(Sig_Gene_Analysis_ID$Sampling_Point)

#Adjusting column names
colnames(Sig_Gene_Analysis_ID)[2] <- "Survival"
colnames(Sig_Gene_Analysis_ID)[4] <- "Sampling Point"

upset4_analysis <- 
  ComplexUpset::upset(Sig_Gene_Analysis_ID, 
                      colnames(Sig_Gene_Analysis_ID)[2:4],
                      sort_intersections_by=c("degree","cardinality"),
                      name = "",
                      width_ratio = 0.1,
                      stripes = c("grey88", "white"),
                      themes = upset_default_themes(text = element_text(size = 15)),
                      matrix = (
                        intersection_matrix(geom = geom_point(shape="circle filled", size =3))
                        + scale_color_manual(
                          values = c("Sampling Point" = "#009E73", "Survival" = "#56B4E9", "Day" = "#D55E00"),
                          guide=guide_legend(override.aes = list(shape = "circle"))
                        )
                      ),
                      queries = list(
                        upset_query(set="Sampling Point", fill="#009E73"),
                        upset_query(set="Survival", fill="#56B4E9"),
                        upset_query(set="Day", fill="#D55E00")
                      )
  )

#Change y-axis labels for teh two bar plots  
upset4_analysis[[2]][["labels"]][["y"]] <- "Number of Genes"
upset4_analysis[[3]][["labels"]][["y"]] <- "Number of Genes"

#Modify the total number of genes in analysis plot (bottom left bar plot)
graph_pract <- upset4_analysis[[3]]
upset4_analysis[[3]] <- 
  graph_pract+
  geom_text(stat = "count", label = c("5,264","5,449","4,449"), nudge_y = 2000, size=4)+ 
  theme(axis.text = element_blank(), axis.title = element_text(size = 8))

upset4_analysis

ggsave("Plots/Fig2_UpsetAnalysis.png", upset4_analysis, width = 7.5, height =7, units = "in")

#
#Clustering Identified Genes using Expression Trajectories####
#Packages:
library(factoextra)
library(tidyverse)
library(rms)
library(RColorBrewer)

#Code expects a certain directory structure with these specific files:
# Master Directory/
#   deseq_results/
#       day/
#         day_Padj_Z_scores.csv
#       sampling_point/
#         sampling_point_Padj_Z_scores.csv
#       survivorship/
#         survivorship_Padj_Z_scores.csv
#       All Analyses/
#         Sig_Genes_AnalysisID.csv

day_z <- read.csv("deseq_results/day/day_Padj_Z_scores.csv", header = TRUE)
survivorship_z <- read.csv("deseq_results/survivorship/survivorship_Padj_Z_scores.csv", header = TRUE)
sampling_point_z <- read.csv("deseq_results/sampling_point/sampling_point_Padj_Z_scores.csv", header = TRUE)
all_sig_genes <- read.csv("deseq_results/All Analyses/Sig_Genes_AnalysisID.csv", header = TRUE)
all_sig_genes <- all_sig_genes[-c(1)]

#Combine all Z scores
z_score <- rbind(day_z, survivorship_z, sampling_point_z) #227430
z_score <- z_score[-c(1)]
#Z scores for shared genes will be the same no matter which analysis so only have
#1 set of z scores for each gene. 
z_score <- unique(z_score) #92130 when divided by 15 = 6142
#Add the analysis identification information
z_score <- merge(z_score, all_sig_genes[,c(1,5)], by=c("gene_id"))


all_clustering_func <- function(Z ,METHOD = "pearson" ,K,ANALYSIS = "All Analyses"){
  Z_scores <- Z
  #Pivot Z scores each gene id is a column, and days are the rows
  Z_scores_wide <- Z_scores %>% pivot_wider(id_cols = c(day), names_from = gene_id,
                                            values_from = Z_score) %>% arrange(day)
  Z_scores_wide <- column_to_rownames(Z_scores_wide, var = "day")
  #Create a table to save which cluster each gene belongs too 
  cluster.correlation <- matrix(0, nrow = ncol(Z_scores_wide), ncol = 1)
  colnames(cluster.correlation) <- c("Cluster")
  rownames(cluster.correlation) <- colnames(Z_scores_wide)
  print("Clustering")
  #Create dissimilarty matrix using pearson correlation.
  cor_diss <- get_dist(t(Z_scores_wide), method = METHOD)
  #Do hierchical clustering with the matrix hclust(cor_diss), and then cut it
  #into how many groups as designated by k. 
  cluster.correlation[,1] <- cutree(hclust(cor_diss), k=K)
  cluster.correlation <- as.data.frame(cluster.correlation)
  cluster.correlation <- rownames_to_column(cluster.correlation, "gene_id")
  
  #Merge with the Z scores the cluster information and save
  Z_scores_grouped <- merge(Z_scores, cluster.correlation,by="gene_id")

  #Plot our clusters
  print("dendrogram plot")
  hclustering <- hclust(cor_diss)
  Z_scores_grouped$Cluster <- factor(Z_scores_grouped$Cluster)
  pdf_name <- paste("deseq_results/All Analyses/Combo_Analysis_clustering_k",K,"_plot.pdf",sep = "")
  pdf(pdf_name, width = 10, height = 10)
  #Created a dendrogram of our clusters
  #Note this takes a very long time to create ~18 minutes
  cluster_dendrogram <- fviz_dend(hclustering, show_labels = FALSE, k=K,palette = "jco",
                                  main = paste("All Analyses clustering k = ", K, sep=""))
  dendrogram_name <- paste("deseq_results/All Analyses/Combo_Analysis_dendrogram_k",K,".rds",sep = "")
  #Saving the dendrogram so it can be used/edited in the future
  saveRDS(cluster_dendrogram, dendrogram_name)
  #Get the default colors used in the dendrogram
  cluster_cols <- cluster_dendrogram[["plot_env"]][["data"]][["labels"]]
  gene_col <- cluster_cols[,3:4]
  colnames(gene_col)[1] <- c("gene_id")
  Z_scores_grouped <- merge(Z_scores_grouped, gene_col, by=c("gene_id"))
  #Cluster num given in cuttree and the order in dendrogram do not match. Using the colors
  #we can get the order of the clusters as seen in the dendrogram left to right
  color_order <- cbind(unique(gene_col$col),1:28)
  color_order <- as.data.frame(color_order)
  colnames(color_order) <- c("col", "Dendo_clus")
  Z_scores_grouped <- merge(Z_scores_grouped, color_order, by=c("col"))
  Z_scores_grouped$Dendo_clus <- as.numeric(Z_scores_grouped$Dendo_clus)
  Z_scores_grouped_filename <- paste("deseq_results/All Analyses/Combo_Analysis_Clustered_Z_scores_k",K,".csv", sep = "")
  write.csv(Z_scores_grouped, Z_scores_grouped_filename)
  colors <- c(brewer.pal(7,"Set2"))
  analysis.id.colors <- c(SDP = colors[1], P = colors[2], SD = colors[3],DP = colors[4],D =colors[5],
                          S = colors[6], SP = colors[7])
  print(cluster_dendrogram)
  #Plot each cluster's gene's expression trajectroy as well as rep curve.
  #Indv genes are colored based off their analysis id. 
  #Clusters are in order as seen in dendrogram
  for (c in 1:K) {
    plot_title <- paste("All Analyses clustering C:",c, sep = "")
    plot_data <- Z_scores_grouped %>% filter(Dendo_clus == c)
    c_color <- unique(plot_data$col)
    cluster_plot <- plot_data %>% 
      ggplot(aes(x=day, y=Z_score))+
      ylim(-3.5,4)+
      geom_point(aes(color=Analysis_ID),show.legend = TRUE)+
      geom_smooth(method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE, show.legend = FALSE,aes(group=gene_id,color=Analysis_ID))+
      geom_smooth(method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE, show.legend = FALSE, color="black")+
      scale_color_manual(values = analysis.id.colors)+
      ggtitle(plot_title)+
      theme_classic()+
      theme(plot.title = element_text(colour = c_color))
    print(cluster_plot)
  }
  graphics.off()
}
all_clustering_func(z_score, K=28)
#This function will take awhile to run. Creating the dendrogram
#alone takes about 18 minutes


#Cluster Number Determination####
#How we decided to divide our genes into 28 clusters

library(factoextra)
library(tidyverse)
library(rms)

#Code expects a certain directory structure with these specific files:
# Master Directory/
#   deseq_results/
#       day/
#         day_Padj_Z_scores.csv
#       sampling_point/
#         sampling_point_Padj_Z_scores.csv
#       survivorship/
#         survivorship_Padj_Z_scores.csv

#Read in Z scores for significant genes
day_z <- read.csv("deseq_results/day/day_Padj_Z_scores.csv", header = TRUE)
survivorship_z <- read.csv("deseq_results/survivorship/survivorship_Padj_Z_scores.csv", header = TRUE)
sampling_point_z <- read.csv("deseq_results/sampling_point/sampling_point_Padj_Z_scores.csv", header = TRUE)

#Combine all Z scores
z_score <- rbind(day_z, survivorship_z, sampling_point_z) #227430
z_score <- z_score[-c(1)]
#Z scores for shared genes will be the same no matter which analysis so only have
#1 set of z scores for each gene. 
z_score <- unique(z_score) #92130 when divided by 15 = 6142

#Pivot Z scores each gene id is a column, and days are the rows
Z_scores_wide <- z_score %>% pivot_wider(id_cols = c(day), names_from = gene_id,
                                          values_from = Z_score) %>% arrange(day)
Z_scores_wide <- column_to_rownames(Z_scores_wide, "day")

#Cluster our genes with varying numbers of clusters
##Sup Table D####
#Aiming to get as many genes in clusters that have 50-500 genes
gene_count_df <- NULL
cluster_df <- NULL
prev_cluster <- NULL
for (K in 15:40) {
  #Create a table to save which cluster each gene belongs too 
  cluster.correlation <- matrix(0, nrow = ncol(Z_scores_wide), ncol = 1)
  colnames(cluster.correlation) <- c("Cluster")
  rownames(cluster.correlation) <- colnames(Z_scores_wide)
  print("Clustering")
  #Create dissimilarty matrix using pearson correlation.
  cor_diss <- get_dist(t(Z_scores_wide), method = "pearson")
  #Do hierchical clustering with the matrix hclust(cor_diss), and then cut it
  #into how many groups as designated by k. 
  cluster.correlation[,1] <- cutree(hclust(cor_diss), k=K)
  cluster.correlation <- as.data.frame(cluster.correlation)
  cluster.correlation <- rownames_to_column(cluster.correlation, "gene_id")
  #count how many genes are within each cluster and recrod how many total clusters
  gene_count <- cluster.correlation %>% dplyr::count(Cluster)
  gene_count$K_value <- K
  gene_count_df <- rbind(gene_count_df, gene_count)
  #Get the 1st, 2nd and 3rd Quantiles
  quarts <- quantile(gene_count$n, prob=c(.25,.5,.75), type = 7)
  #Put cluster in high, low or right group based on gene size
  hi <- gene_count %>% filter(n > 500)
  lo <- gene_count %>% filter(n < 50)
  ri <- gene_count %>% filter(n < 500 & n >50)
  #For each gene designated wether it is in a high, low, or right cluster
  cluster.correlation$Des <- NA
  cluster.correlation$Des[cluster.correlation$Cluster %in% hi$Cluster] <- "H"
  cluster.correlation$Des[cluster.correlation$Cluster %in% lo$Cluster] <- "L"
  cluster.correlation$Des[cluster.correlation$Cluster %in% ri$Cluster] <- "R"
  
  #Create a table that has summary of cluster sizes, number of cluster in our ideal zone
  #and number of genes in our ideal zone in a row for each total cluster number variation
  Des_count <- cluster.correlation %>% dplyr::count(Des) %>% column_to_rownames(var = "Des")
  if (is.null(prev_cluster)){
    row <- cbind(K, min(gene_count$n),quarts[1],mean(gene_count$n), quarts[2],quarts[3], max(gene_count$n),
                 nrow(lo), nrow(ri), nrow(hi), Des_count["L",], Des_count["R",], Des_count["H",], NA, NA)
    colnames(row) <- c("Cluster", "Min","Q1","Mean","Median","Q3","Max",
                       "Low_Cluster", "In_Range_Cluster", "High_Cluster", "Low_Genes", "In_Range_Genes", "High_Genes",
                       "Genes_Gained", "Genes_Lost")
    cluster_df <- as.data.frame(rbind(cluster_df, row))
    prev_cluster <- cluster.correlation
  }else{
    #When past the first total cluster number, also recording how many genes were gained and lost in our ideal zone
    #compared to previous cluster number
    loprev <- prev_cluster %>% filter(Des == "L")
    lostgene <- cluster.correlation %>% filter(Des == "L" & !gene_id %in% loprev$gene_id)
    riprev <- prev_cluster %>% filter(Des == "R")
    gaingene <- cluster.correlation %>% filter(Des == "R" & !gene_id %in% riprev$gene_id)
    row <- cbind(K, min(gene_count$n),quarts[1],mean(gene_count$n), quarts[2],quarts[3], max(gene_count$n),
                 nrow(lo), nrow(ri), nrow(hi), Des_count["L",], Des_count["R",], Des_count["H",], nrow(gaingene), nrow(lostgene))
    colnames(row) <- c("Cluster", "Min","Q1","Mean","Median","Q3","Max",
                       "Low_Cluster", "In_Range_Cluster", "High_Cluster", "Low_Genes", "In_Range_Genes", "High_Genes",
                       "Genes_Gained", "Genes_Lost")
    cluster_df <- as.data.frame(rbind(cluster_df, row))
    prev_cluster <- cluster.correlation
  }
  rownames(cluster_df) <- 1:nrow(cluster_df)
}

dir.create("deseq_results/Cluster Number Determination")
#Save file and this our Supplementary table C
write.csv(cluster_df, "deseq_results/Cluster Number Determination/K15_40_cluster_df.csv")
#When looking at this table we can see the last massive gain of genes in our ideal zone
#happens when we increase to 27 clusters and we begin to loose genes when we increase
#to cluster 30. This put our ideal cluster number between 27 and 29. To determine
#which number we should use we looked at the size of the clusters being broken apart between
#these transition

gene_count_df %>% filter(K_value == 27 | K_value == 28 | K_value == 29)
#Clusters being broken were just manually identified by comparing cluster numbers
#between different K values

#Visualize the changing clusters
#Mark clusters that are going being broken apart with red, and the ones that 
#are the new clusters in green
##Sup Fig 2  Cluster Determination####
gene_count_df$color <- "NC"
gene_count_df$color[gene_count_df$K_value == 27 & gene_count_df$n == 458] <- "R"
gene_count_df$color[gene_count_df$K_value == 28 & gene_count_df$n == 167] <- "B"
gene_count_df$color[gene_count_df$K_value == 28 & gene_count_df$n == 130] <- "R"
gene_count_df$color[gene_count_df$K_value == 28 & gene_count_df$n == 328] <- "R"
gene_count_df$color[gene_count_df$K_value == 29 & gene_count_df$n == 66] <- "B"
gene_count_df$color[gene_count_df$K_value == 29 & gene_count_df$n == 101] <- "B"

gene_count_df$point <- "No Change"
gene_count_df$point[gene_count_df$K_value == 27 & gene_count_df$n == 458] <- "Splitting Cluster"
gene_count_df$point[gene_count_df$K_value == 28 & gene_count_df$n == 167] <- "Splitting Cluster"
gene_count_df$point[gene_count_df$K_value == 28 & gene_count_df$n == 130] <- "New Cluster"
gene_count_df$point[gene_count_df$K_value == 28 & gene_count_df$n == 328] <- "New Cluster"
gene_count_df$point[gene_count_df$K_value == 29 & gene_count_df$n == 66] <- "New Cluster"
gene_count_df$point[gene_count_df$K_value == 29 & gene_count_df$n == 101] <- "New Cluster"

#2 clusters which have 66 genes
gene_count_df[324,4] <- "NC"
gene_count_df[324,5] <- "No Change"

gene_count_df$'Cluster Change:' <- gene_count_df$point

gene_distribute_graph<- 
  gene_count_df %>%filter(K_value == "27" |K_value == "28" |K_value == "29") %>% 
  ggplot(aes(x=K_value, y=n, group=K_value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color = color, shape = `Cluster Change:`, size=`Cluster Change:`, alpha=`Cluster Change:`), fill="white", position = position_jitter(0.075))+ #size=2
  coord_flip()+
  scale_shape_manual(values = c(21,20,19))+
  scale_color_manual(values = c("dodgerblue", "black", "red"), guide="none")+
  scale_size_manual(values = c(2.5,2,2.5))+
  scale_alpha_manual(values = c(1, 0.5, 1))+
  ylab("Number of Genes")+
  xlab("Total Number of Clusters")+
  scale_x_reverse(breaks=c(27,28,29))+
  ggtitle("Number of Genes per Cluster")+
  theme_bw()+
  theme(text = element_text(size=10), legend.position = "bottom")

ggsave("deseq_results/Cluster Number Determination/SupFig2_Gene_Distribution.png", gene_distribute_graph, height = 6, width = 7, units = "in")

#
#Classifying Cluster Trajectories####
#Want to determine if a cluster's trajectory is linear or complex. 
library(tidyverse)
library(Rmisc)
library(rms)

#Code expects a certain directory structure and files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       Combo_Analysis_Clustered_Z_scores_k28.csv

#Will need first to get the z scores for our representative curves
z_scores <- read.csv("deseq_results/All Analyses/Combo_Analysis_Clustered_Z_scores_k28.csv", header = TRUE)
z_scores <- z_scores[-c(1)]

#For every cluster we are getting the mean z score for each day across all genes in that cluster
rep_z <- NULL
for (c in unique(z_scores$Cluster)) {
  clus_z <- z_scores %>% filter(Cluster == c)
  clus_z_summary <- summarySE(data = clus_z, measurevar = "Z_score", groupvars = c("day","survivorship","Cluster","Dendo_clus"))
  rep_z <- rbind(rep_z, clus_z_summary)
}
dim(rep_z)
# [1] 420   9
# 420/15 28

#To determine if a cluster is complex or linear we are going to test if a curve is linear
#by doing a correlation test between the z scores of the rep line and time.     

##Sup Table E####
#Add in sampling point column
rep_z$SP <- rep(1:15, 28)

lm_results <- NULL
for (c in 1:28) {
  #Run three seperate linear regression for each of the analyses
  day <- lm(Z_score~day, data = rep_z %>% filter(Dendo_clus == c))
  surv <- lm(Z_score~survivorship, data = rep_z %>% filter(Dendo_clus == c))
  sp <- lm(Z_score~SP, data = rep_z %>% filter(Dendo_clus == c))
  aov_day <- anova(day)
  aov_surv <- anova(surv)
  aov_sp <- anova(sp)
  lm_results <- rbind(lm_results,c(c,aov_day$`Pr(>F)`[1], day$coefficients[2]),
                      c(c,aov_surv$`Pr(>F)`[1], surv$coefficients[2]*-1),
                      c(c,aov_sp$`Pr(>F)`[1], sp$coefficients[2])
  )
}
lm_results <- as.data.frame(lm_results)
lm_results$Analysis <- rep(c("Day","Survival","Sampling Point"),28)
colnames(lm_results)[1:3] <- c("Dendorogram_Cluster_Number", "Pvalue", "Coefficient")
#Designate Clusters as Complex, LinearUp or LinearDown
#First designate all clusters as complex
lm_results$Designation <- "Complex"
#If a cluster has a p-value less than 0.05/28 it gets designated as Linear. Its
#direction of up or down is then determined the sign of its coefficient. 
lm_results$Designation[lm_results$Pvalue < 0.05/28 & lm_results$Coefficient > 0] <- "LinearUp"
lm_results$Designation[lm_results$Pvalue < 0.05/28 & lm_results$Coefficient < 0] <- "LinearDown"

#Create Cluster's Official ID as its Designation and what iteration of that designation
#it is. 
c <- 1
lu <- 1
ld <- 1
for (n in 1:28) {
  dendo_desg <- unique(lm_results %>% filter(Dendorogram_Cluster_Number == n) %>% pull(Designation))
  if(dendo_desg == "Complex"){
    official_id <- paste(dendo_desg,"-",c,sep = "")
    c = c+1
  }else if(dendo_desg == "LinearUp"){
    official_id <- paste(dendo_desg,"-",lu,sep = "")
    lu = lu+1
  }else if(dendo_desg == "LinearDown"){
    official_id <- paste(dendo_desg,"-",ld,sep = "")
    ld = ld+1
  }
  lm_results$Official_ID[lm_results$Dendorogram_Cluster_Number == n] <- official_id
}

#Save Table. This is Supplementary Table E
write.csv(lm_results,"deseq_results/All Analyses/Cluster_Designation_Linear_Regression.csv")

color_clus_id <- read.csv("deseq_results/All Analyses/cluster_color_id_ref.csv", header = TRUE)
lm_results <- merge(lm_results, color_clus_id[c(3,5)], by=c("Dendo_clus"))






#Additional Figures####
#Figures in the paper that use the data generated in the previous sections, but
#are not directly made in those sections

##Figure 3 Simplified Dendrogram####
#Want to simplify dendrogram so each cluster is represented by one colored line,
#labeled at the end, and colors are more distinct then default colors used in
#Clustering Identified Genes section
library(RColorBrewer)
library(tidyverse)
library(factoextra)

#Code expects a certain directory structure and files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       Combo_Analysis_Clustered_Z_scores_k28.csv
#       Combo_Analysis_dendrogram_k28.rds

dendo_28 <- readRDS("deseq_results/All Analyses/Combo_Analysis_dendrogram_k28.rds")
col_info <- read.csv("deseq_results/All Analyses/Combo_Analysis_Clustered_Z_scores_k28.csv", header = TRUE)
dendo_white <- dendo_28
colors <- unique(col_info$col)
colors <- append(colors, "black")

#We will being pull out sections of the dendogram's data, edit it, and then 
#replace the original data with the edited version. 

#The dendogram is basically a series of vertical and horizontal lines plotted onto
#a graph. This data tell tell us the x and y coordinates for each lines starting and end
#point, as well as its colors, linewidth and linetype.
df <- dendo_white[["layers"]][[1]][["data"]]

#Each color is represented by an unique color and our goal is to only 
#have one vertical line for each cluster

df_recolored <- NULL
for (c in colors) {
  if(c == "black"){ #Keep all lines that are black
    df_recolored <- rbind(df_recolored, df %>% filter(col ==c))
  }else{
    #Pull out lines with that colors
    df_color <- df %>% filter(col == c)
    #Take the y value in the first row. This is the first horizontal for that
    #cluster
    color_y <- df_color[1,2]
    #filter rows so only those whose y value is greater or equal color_y are kept.
    #This will only keep two rows. The first horizontal line row and the first vertical.
    df_color <- df_color %>% filter(y >= color_y)
    #Change the y end value for the vertical line, so all lines end at the same place
    df_color[2,4] <- 1
    df_recolored <- rbind(df_recolored, df_color)
  }
}

#Replace the original data with the edited data
dendo_white[["layers"]][[1]][["data"]] <- df_recolored
#Should see only 1 vertical line for each cluster and they all ending at 1. 
dendo_white

#Adding labels to dendrogram
#Going take the labels from an example graph and add them into our graph
df <- scale(USArrests)
res.hc <- hclust(dist(df))
dend <- fviz_dend(res.hc)
test_label <- dend[["layers"]][[2]]
dendo_white[["layers"]][[3]] <- test_label

#Create our labels
deno_labels <- factor(c("Complex-1","Complex-2","Complex-3","Complex-4","Complex-5","Complex-6",
                        "LinearUp-1",
                        "Complex-7",
                        "LinearUp-2","LinearUp-3","LinearUp-4", "LinearUp-5",
                        "Complex-8",
                        "LinearUp-6", "LinearUp-7","LinearUp-8",
                        "Complex-9", "Complex-10","Complex-11",
                        "LinearDown-1", "LinearDown-2", "LinearDown-3", "LinearDown-4","LinearDown-5",
                        "Complex-12",
                        "LinearDown-6",
                        "Complex-13","Complex-14"), levels = c("Complex-1","Complex-2","Complex-3","Complex-4","Complex-5","Complex-6",
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

#Pull our temporary labels, keep first 28 rows, and then replace labels with our
#labels
dendo_label_edit <- dendo_white[["layers"]][[3]][["data"]]
dendo_label_edit <- dendo_label_edit[1:28,]
dendo_label_edit$label <- deno_labels

#Now need to get the x position for our labels 
#For each color that isn't black, take the row that is a vertical line and
#pull out the x_position of that line
x_pos <- NULL
for (c in colors) {
  if(c != "black"){
    df_color <- dendo_white[["layers"]][[1]][["data"]] %>% filter(col == c)
    for (r in 1:nrow(df_color)) {
      row <- df_color[r,]
      if(row$x == row$xend){
        x_pos <- append(x_pos, row$x)
      }
    }
  }
}
#Sort so its in increasing order
x_pos <- sort(x_pos)
#Replace the x values with new x positions
dendo_label_edit$x <- x_pos
dendo_label_edit
dendo_white[["layers"]][[3]][["data"]] <- dendo_label_edit
dendo_white
#Now have labels at the bottom of the graph. The first 6 rows are on top of each
#other and not readable. We will address this later on. 

#Recolor the clusters
#Want to have distinct colors for our clusters. Not concerned about each cluster
#having a unique color, just so that it easy to tell them apart from each other.
display.brewer.all(type="qual", colorblindFriendly = TRUE)
#Pull out the first 7 colors fro Dark 2 pallet and then repeat it 4 times so we
#have 28 color values
colors <- brewer.pal(7,"Dark2")
colors <- rep(colors, 4)
#Associate the colors with a cluster
col_info <- unique(col_info[c(2,15:16)])
col_info <- col_info %>% arrange(Dendo_clus)
col_info$New_colors <- colors
#Replace the og colors with the new colors
for (c in col_info$col) {
  df_recolored$col[df_recolored$col == c] <- col_info %>% filter(col == c) %>% pull(New_colors)
}
dendo_white[["layers"]][[1]][["data"]] <- df_recolored
dendo_white

#Reorientation and label adjustments
dendo_label_edit$y <- 0.95
dendo_label_edit$angle <- 0
dendo_label_edit$angle <- 0
dendo_label_edit$hjust <- 0

#For the first 5 clusters, instead of each having a label we are going create
#a label to represent all 5. 
#Get the average x position of first 5 clusters
x_avg <- mean(dendo_label_edit[1:5,1])
#Create new row whose x psotion in the average x position
dendo_label_edit[29,] <- dendo_label_edit[1,]
dendo_label_edit[29,1] <- x_avg
#Remove the first 5 individual rows
dendo_label_edit <- dendo_label_edit[6:29,]
#Change the new rows label so it represents all 5 of the clusters
dendo_label_edit$label <- as.character(dendo_label_edit$label)
dendo_label_edit[24,3] <- c("Complex-1:5")

#Further adjustments so labels do not overlap
dendo_label_edit$cex <- 3
dendo_label_edit[23,1] <- 6100
dendo_label_edit[12,1] <- 2958
dendo_label_edit[13,1] <- 3049
dendo_label_edit[8,1] <- 2300

dendo_white[["layers"]][[3]][["data"]] <- dendo_label_edit
dendo_white+
  coord_flip(ylim = c(2,0.75))+
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        title = element_blank())

ggsave("Plots/Fig3_Dendo_Vertical.png", dendo_white+
         coord_flip(ylim = c(2,0.7))+
         theme(axis.ticks = element_blank(),
               axis.title = element_blank(),
               axis.text.y = element_blank(),
               title = element_blank()), width = 4, height = 8, units = "in")

###Color ID Sheet####
#Create reference sheet with the names and colors for the clusters and genes
library(tidyverse)
library(RColorBrewer)

#Code expects a certain directory structure and files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       Combo_Analysis_Clustered_Z_scores_k28.csv

color_id <- read.csv("deseq_results/All Analyses/Combo_Analysis_Clustered_Z_scores_k28.csv", header = TRUE)
gene_color_id <- unique(color_id[,c(3,14:16,2)])
cluster_color_id <- unique(color_id[,c(15:16,2)])
cluster_color_id <- cluster_color_id %>% arrange(Dendo_clus)
cluster_color_id$Official_ID <- c("Complex-1","Complex-2","Complex-3","Complex-4","Complex-5","Complex-6",
                                  "LinearUp-1",
                                  "Complex-7",
                                  "LinearUp-2","LinearUp-3","LinearUp-4", "LinearUp-5",
                                  "Complex-8",
                                  "LinearUp-6", "LinearUp-7","LinearUp-8",
                                  "Complex-9", "Complex-10","Complex-11",
                                  "LinearDown-1", "LinearDown-2", "LinearDown-3", "LinearDown-4","LinearDown-5",
                                  "Complex-12",
                                  "LinearDown-6",
                                  "Complex-13","Complex-14")

display.brewer.all(type="qual", colorblindFriendly = TRUE)
colors <- brewer.pal(7,"Dark2")
cluster_color_id$New_color <- colors
gene_color_id <- merge(gene_color_id, cluster_color_id[c(1,4,5)], by="Cluster")
write.csv(cluster_color_id,"deseq_results/All Analyses/cluster_color_id_ref.csv")
write.csv(gene_color_id,"deseq_results/All Analyses/gene_color_id_ref.csv")

#
##Sup Figure 3 Recolored Dendrogram####
#Complete dendrogram with each gene plotted but recolored and labeled to match 
#Simplified dendrogram 
library(RColorBrewer)
library(tidyverse)
library(factoextra)

#Code expects a certain directory structure and files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       Combo_Analysis_dendrogram_k28.rds
#       cluster_color_id_ref.csv

dendo_28 <- readRDS("deseq_results/All Analyses/Combo_Analysis_dendrogram_k28.rds")
dendo_supp <- dendo_28
cluster_col_info <- read.csv("deseq_results/All Analyses/cluster_color_id_ref.csv", header = TRUE)

#Create labels
df <- scale(USArrests)
res.hc <- hclust(dist(df))
dend <- fviz_dend(res.hc)
test_label <- dend[["layers"]][[2]]
dendo_supp[["layers"]][[3]] <- test_label
deno_labels <- factor(c("Complex-1","Complex-2","Complex-3","Complex-4","Complex-5","Complex-6",
                        "LinearUp-1",
                        "Complex-7",
                        "LinearUp-2","LinearUp-3","LinearUp-4", "LinearUp-5",
                        "Complex-8",
                        "LinearUp-6", "LinearUp-7","LinearUp-8",
                        "Complex-9", "Complex-10","Complex-11",
                        "LinearDown-1", "LinearDown-2", "LinearDown-3", "LinearDown-4","LinearDown-5",
                        "Complex-12",
                        "LinearDown-6",
                        "Complex-13","Complex-14"), levels = c("Complex-1","Complex-2","Complex-3","Complex-4","Complex-5","Complex-6",
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
dendo_label_edit <- dendo_supp[["layers"]][[3]][["data"]]
dendo_label_edit <- dendo_label_edit[1:28,]
dendo_label_edit$label <- deno_labels

#Get the X_position for the labels
#Similar to simplified version's, except now we are getting the mean
#x position of all the vertical lines in a cluster
x_pos <- NULL
for (c in unique(col_info$col)) {
  if(c != "black"){
    df_color <- dendo_supp[["layers"]][[1]][["data"]] %>% filter(col == c)
    color_x_pos <- NULL
    for (r in 1:nrow(df_color)) {
      row <- df_color[r,]
      if(row$x == row$xend){
        color_x_pos <- append(color_x_pos, row$x)
      }
    }
    x_pos <- append(x_pos, mean(color_x_pos))
  }
}
dendo_label_edit$x <- x_pos
dendo_label_edit$y <- -0.025

#Create the rep label for first 5 clusters
dendo_label_edit[29,] <- dendo_label_edit[1,]
dendo_label_edit[29,1] <- mean(dendo_label_edit[1:5,1])
dendo_label_edit$label <- as.character(dendo_label_edit$label)
dendo_label_edit[29,3] <- "Complex-1:5"
dendo_label_edit <- dendo_label_edit[6:29,]
dendo_supp[["layers"]][[3]][["data"]] <- dendo_label_edit

#Change cluster colors
df <- dendo_supp[["layers"]][[1]][["data"]]
for (c in unique(cluster_col_info$col)) {
  df$col[df$col == c]<- cluster_col_info %>% filter(col == c) %>% pull(New_color)
}
dendo_supp[["layers"]][[1]][["data"]] <- df
dendo_supp+theme(title = element_blank())

ggsave("Plots/SupFig3_Dendo.png", dendo_supp+theme(title = element_blank()), width = 11, height = 8, units = "in") 

##Figure 4 Rep Curves####
library(tidyverse)
library(rms)
library(RColorBrewer)

#Code expects a certain directory structure and files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       Combo_Analysis_Clustered_Z_scores_k28.csv
#       Cluster_Designation_Linear_Regression.csv

#Read in Z scores for all genes, and for rep curves which has the cluster's trajectory
#Designation
rep_z <- read.csv("deseq_results/All Analyses/Cluster_Designation_Linear_Regression.csv", header = TRUE)
rep_z <- rep_z[-c(1)]
z_scores <- read.csv("deseq_results/All Analyses/Combo_Analysis_Clustered_Z_scores_k28.csv", header = TRUE)
z_scores <- z_scores[-c(1)]

#Clusters are named for their trajectory designation (Complex, LinearDown, LinearUp) and then
#what iteration they are of that designation going from left to right or bottom to top
#on the dendrogram depending on its orientation. 
cluster_names <- unique(rep_z[c(1,6)])
colnames(cluster_names)[1] <- "Dendo_clus"
z_scores_names <- merge(z_scores, cluster_names, by=c("Dendo_clus"))

#Factor the names so complex trajectories are plotted first, than LinearUp and then LinearDown
z_scores_names$Official_ID <- factor(z_scores_names$Official_ID, levels = c("Complex-1","Complex-2","Complex-3","Complex-4",
                                                              "Complex-5","Complex-6","Complex-7","Complex-8",
                                                              "Complex-9", "Complex-10", "Complex-11","Complex-12",
                                                              "Complex-13","Complex-14", "LinearUp-1","LinearUp-2",
                                                              "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                                              "LinearUp-7", "LinearUp-8", "LinearDown-1", "LinearDown-2",
                                                              "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6"))

z_scores_names <- z_scores_names %>% arrange(Official_ID)
#Create Labels for the clusters which includes the cluster names and the number of genes
#in that cluster
cluster_label <- NULL
for (c in unique(z_scores_names$Official_ID)) {
  print(c)
  num_genes <- nrow(unique(z_scores_names[,c(3,16)] %>% filter(Official_ID == c)))
  cluster_label <- append(cluster_label,
                          paste(c,"\n",num_genes," genes",sep = ""))
}
#Refactorize names, but now include the labels
z_scores_names$Official_ID <- factor(z_scores_names$Official_ID, levels = c("Complex-1","Complex-2","Complex-3","Complex-4",
                                                              "Complex-5","Complex-6","Complex-7","Complex-8",
                                                              "Complex-9", "Complex-10", "Complex-11","Complex-12",
                                                              "Complex-13","Complex-14", "LinearUp-1","LinearUp-2",
                                                              "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                                              "LinearUp-7", "LinearUp-8", "LinearDown-1", "LinearDown-2",
                                                              "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6"),
                              labels = cluster_label)

#Want to color the curves in the same colors as the clusters are in the dendrogram
display.brewer.all(type="qual", colorblindFriendly = TRUE)
colors <- brewer.pal(7,"Dark2")
colors <- rep(colors, 4)
colors <- as.data.frame(cbind(1:28,colors))
colnames(colors)[1] <- c("Dendo_clus")

colors <- merge(unique(z_scores_names[,c(1,16)]), colors, by=c("Dendo_clus"))
colors <- colors %>% arrange(Official_ID)

z_scores_names$color[gsub("-.*","",z_scores_names$Official_ID) == "Complex"] <- "white"
z_scores_names$color[gsub("-.*","",z_scores_names$Official_ID) == "LinearDown"] <- "gray65"
z_scores_names$color[gsub("-.*","",z_scores_names$Official_ID) == "LinearUp"] <- "gray90"


#Plot the rep curves
rep_curves <- 
  z_scores_names %>% ggplot(aes(x=day, y=Z_score))+
  ylab("Z Score")+
  xlab("Day")+
  facet_wrap(~Official_ID, nrow=4)+
  geom_smooth(method = "lm",aes(color=Official_ID),formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))),show.legend = FALSE)+
  scale_color_manual(values = colors$colors)+
  theme_bw()+
  coord_cartesian(clip="off", ylim=c(-2, 3)) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=3.2, ymax=4.75),alpha=0.4,fill=z_scores_names$color)+
  theme(text = element_text(size=7),
        strip.background = element_rect(fill=NA))
  
ggsave("Plots/Fig4_RepCurves.png",rep_curves, width = 6, height = 6, units = "in")
