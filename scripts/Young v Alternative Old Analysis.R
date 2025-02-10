#Introduction####
#Title: Alternative YvO analyses comparison

#Purpose: Redo the YvO analysis with alternative old time points.
#Using two sequential timepoints for our old time points (exception being
#Days 3 and 6 which will still stand for our young timepoint). We want
#to see how the changing the old cutoff point impact what genes are identified. 
#Will show how to create Supplementary Figures 8 and 9 

#Setup:
#Each section begins with brief description of the goal, what R packages are 
#required and the directory structure including files the code expects. Assumes 
#you have made gene_color_id_ref.csv (See Identifying_and_Clustering.R) and
#gse_salmon_tximeta_yvo_AltOld#v36.RDS

#Created: 8/19/2024, KMH
#Last Edited: 11/08/24, KMH

#Deseq Gene Identification####
library(BiocManager)
library(DESeq2)
library(tidyverse)
library(Rmisc)

#Code expects a certain directory structure with this specific files:
# Master Directory/
#   Aging Analysis Results/
#     Old v Young Sim/
#       Alternative old Code/
#         gse_salmon_tximeta_yvo_AltOld#v36.RDS (11x)

#Note: the AltOld# is the number days for the old timepoints. 
#So code that is uding Days 14 and 17 as their old time points would
#be listed as gse_salmon_tximeta_yvo_1417v36.RDS

#Run our deseq analysis for each of the tximeta files which end in RDS
deseq_gene <- function(DIR){
  #Get list of files in designated directory that ends with RDS
  RDS_list <- list.files(path = DIR, pattern = "RDS")
  #For each of the files run our Deseq code
  for (r in RDS_list) {
    gse <- readRDS(paste(DIR, r, sep = ""))
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
    #Created a new directory which is named by AltOld#v36 in the RDS files name
    new_dir <- paste(DIR, substr(r, 24, nchar(r)-4),"/", sep = "")
    dir.create(new_dir)
    #Save files
    write.csv(as.data.frame(sig_test), paste(new_dir, "OvY_PadjSig_deseq_",substr(r, 24, nchar(r)-4),".csv", sep = ""))
    write.csv(as.data.frame(res), paste(new_dir, "OvY_All_deseq_", substr(r, 24, nchar(r)-4),".csv", sep = ""))
    read_counts <- counts(dds, normalized=TRUE)
    write.csv(read_counts, paste(new_dir,"OvY_Norm_counts_deseq_", substr(r, 24, nchar(r)-4), ".csv", sep = ""))
  }
}

dir <- paste("Aging Analysis Results/Old v Young Sim/Alternative old Code/")
deseq_gene(dir)

#Shared Genes with Day59 old analysis and Multi timepoint analysis####
library(tidyverse)

#Code expects a certain directory structure with this specific files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       gene_color_id_ref.csv
#   Aging Analysis Results/
#     Old v Young Sim/
#       OvY_PadjSig_deseq.csv  
#       Alternative old Code/
#         AltOld#v36/
#           OvY_PadjSig_deseq_AltOld#v36.csv

shared_genes <- function(DIR){
  #Read in identified genes from Old59 and multi-timepoint analysis
  OG_ovy <- read.csv("Aging Analysis Results/Old v Young Sim/OvY_PadjSig_deseq.csv", header = TRUE)
  OG_ovy <- OG_ovy[-c(1)]
  multi_tp <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv")
  multi_tp <- multi_tp[-c(1)]
  summary_table <- NULL
  #For every directory in our designated directory
  for (d in list.dirs(DIR)[-1]) {
    #Read in the PadjSig genes for that altnerative old anlaysis
    alt_ovy <- read.csv(paste(d,"/",list.files(path = d, "PadjSig"),sep = ""), header = TRUE)
    #Pull out shared genes from multi-timepoint and Old59 analyses
    shared_genes <- OG_ovy %>% filter(gene_id %in% alt_ovy$gene_id)
    share_genes_multi <- multi_tp %>% filter(gene_id %in% alt_ovy$gene_id)
    id <- substr(d, 61, nchar(d))
    #Create a row for summary tables that infoirmation regarding number of shared genes
    #number of unique genes, and total number of genes
    new_row <- as.data.frame(cbind(id, nrow(shared_genes), nrow(share_genes_multi),
                                   nrow(OG_ovy),nrow(multi_tp),nrow(alt_ovy),
                                   nrow(OG_ovy)-nrow(shared_genes), nrow(alt_ovy)-nrow(shared_genes),
                                   nrow(multi_tp)-nrow(share_genes_multi), nrow(alt_ovy)-nrow(share_genes_multi)))
    colnames(new_row) <- c("Old_ID", "Num_shared_w_Old59", "Num_shared_w_MultiTP",
                           "Total_Old59","Total_MultiTP","Total_Old_Alt", 
                           "Unique_Old59", "Unique_Old_Alt_v_Old59",
                           "Unique_MultiTP", "Unique_Old_Alt_v_MultiTP")
    summary_table <- rbind(summary_table, new_row)
    #Save shared genes table to the alternative old's directory
    write.csv(shared_genes, paste(d,"/shared_genes_w_Old59_",id,".csv", sep = ""))
    write.csv(share_genes_multi, paste(d,"/shared_genes_w_MultiTP_",id,".csv", sep = ""))
  }
  #save the summary table
  write.csv(summary_table, "Aging Analysis Results/Old v Young Sim/Alternative old Code/shared_gene_summary.csv")
}

dir <- paste("Aging Analysis Results/Old v Young Sim/Alternative old Code")
shared_genes(dir)

#YvO Direction Calls####
#Get the direction calls for our alternative old analyses

library(tidyverse)

#Code expects a certain directory structure with this specific files:
# Master Directory/
#   Aging Analysis Results/
#     Old v Young Sim/
#       Alternative old Code/
#         AltOld#v36/
#           OvY_PadjSig_deseq_AltOld#v36.csv
#           OvY_Norm_counts_deseq_AltOld#v36.csv

direction_calc_per_analysis <- function(DIR){
  for (d in list.dirs(DIR)[-1]) {
    id <- substr(d, 61, nchar(d))
    print
    alt_ovy <- read.csv(paste(d,"/",list.files(path = d, "PadjSig"),sep = ""), header = TRUE)
    alt_count <- read.csv(paste(d,"/",list.files(path = d, "Norm"), sep = ""), header = TRUE)
    Sig_YVO_Counts <- alt_count %>% filter(X %in% alt_ovy$gene_id)
    colnames(Sig_YVO_Counts)[1] <- "gene_id"
    Sig_YVO_Counts$Young_Count <- NA
    Sig_YVO_Counts$Old_Count <- NA
    Sig_YVO_Counts$Direction <- NA
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
    Sig_YVO_Counts$Analysis_ID <- id
    write.csv(Sig_YVO_Counts, paste(d,"/Sig_Gene_Direction_Change_", id,".csv", sep = ""))
  }
}

dir <- paste("Aging Analysis Results/Old v Young Sim/Alternative old Code")
direction_calc_per_analysis(dir)

#Comparisons with Day59 Old Analysis####
library(tidyverse)
library(RColorBrewer)

#Code expects a certain directory structure with this specific files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       gene_color_id_ref.csv
#   Aging Analysis Results/
#     Sig_Gene_Direction_Change_59v36.csv
#     Old v Young Sim/
#       Alternative old Code/
#         AltOld#v36/
#           Sig_Gene_Direction_Change_AltOld#v36.csv

O59_Comparison <- function(DIR){
  O59 <- read.csv("Aging Analysis Results/Old v Young Sim/Sig_Gene_Direction_Change_59v36.csv", header = TRUE)
  O59_short <- O59[c(2,17)]
  colnames(O59_short)[2] <- c("O59_Direction")
  master_list <- NULL
  for (d in list.dirs(DIR)[-1]){
    direction_calls <- read.csv(paste(d,"/", list.files(path = d, "Direction"), sep = ""), header = TRUE)
    direction_calls <- direction_calls[c(2,17:18)]
    direction_calls_59 <- full_join(direction_calls, O59_short, by=c("gene_id"))
    id <- substr(d, 61, nchar(d))
    direction_calls_59$Analysis_ID <- id
    master_list <- rbind(master_list, direction_calls_59)
  }
  print(master_list)
}

dir <- paste("Aging Analysis Results/Old v Young Sim/Alternative old Code")

D59_compare <- O59_Comparison(dir)

#Create a comparison column that designates genes
D59_compare$Comparison <- NA
for (r in 1:nrow(D59_compare)) {
  if(is.na(D59_compare[r,2])){
    D59_compare[r,5] <- "Unique_Old59" #Found in Day59 old analysis, not alternative
  }else if (is.na(D59_compare[r,4])){
    D59_compare[r,5] <- "Unique_OldAlt" #Found in Alternative old analysis, not Day59
  }else if (D59_compare[r,2] == D59_compare[r,4]){
    D59_compare[r,5] <- "Direct_Call_Match" #Found in both and their calls match
  }else if (D59_compare[r,2] != D59_compare[r,4]){
    D59_compare[r,5] <- "Direct_Call_Differ" #Found in both and their calls differ
  }
}

##SupFig 8 All identified genes comparison####
#Get the comparison counts for each analysis ID
D59_compare_count <- D59_compare %>% dplyr::count(Analysis_ID, Comparison)
D59_compare_count$Comparison <- factor(D59_compare_count$Comparison, levels = c("Unique_Old59","Unique_OldAlt","Direct_Call_Differ","Direct_Call_Match"))
D59_compare_count <- D59_compare_count %>% arrange(desc(Comparison))
D59_labeled_count <- NULL
for (a in unique(D59_compare_count$Analysis_ID)) {
  a_count <- D59_compare_count %>% filter(Analysis_ID == a)
  for (r in 1:nrow(a_count)) {
    if (r == 1){
      a_count[r,4] <- a_count[r,3]/2
    }else{
      a_count[r,4] <- sum(a_count[1:r-1,3])+(a_count[r,3]/2)
    }
  }
  D59_labeled_count <- rbind(D59_labeled_count, a_count)
}

#label position adjustments
D59_labeled_count[3,4] <- D59_labeled_count[3,4]+100
D59_labeled_count[7,4] <- D59_labeled_count[7,4]+55
D59_labeled_count[11,4] <- D59_labeled_count[11,4]+30
D59_labeled_count[15,4] <- D59_labeled_count[15,4]+30
D59_labeled_count[27,4] <- D59_labeled_count[27,4]+30

D59_all_cnt_plot <- 
  D59_labeled_count %>% ggplot()+
  geom_bar(aes(x=Analysis_ID, y=n, fill=Comparison),stat="identity")+
  geom_text(aes(x=Analysis_ID, y=V4, label=n))+
  scale_fill_manual(values = c(brewer.pal(8,"Paired")[c(1,3,5,7)]),
                    labels = c("Old = Day59 Analysis Unique",
                               "Old = Alternative Analysis Unique",
                               "Shared Genes, Calls Differ",
                               "Shared Genes, Calls Match"))+
  scale_x_discrete(labels = c("D10 & 14", "D14 & 17", "D17 & 23", "D23 & 27",
                              "D27 & 31", "D31 & 36", "D36 & 38", "D38 & 42",
                              "D42 & 48", "D48 & 50", "D50 & 55"))+
  ylab("Number of Genes")+
  xlab("Alternative Old Analyses")+
  guides(fill=guide_legend(nrow = 2, byrow=TRUE))+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 325, hjust = 0.15))

ggsave("Plots/SupFig8_Comparison_D59_total_Genes.png", D59_all_cnt_plot, height = 6, width = 6, units = "in")

##SupFig 9 Trajectory Designation Copare####
Multi_Tp_sig <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)
#Add information about genes clusters
D59_w_multi <- merge(Multi_Tp_sig[c(3,7)], D59_compare, by=c("gene_id"))

#Add on multi-timepoint genes that were identified in the alternative or O59 analyses
for (a in unique(D59_w_multi$Analysis_ID)) {
  multi_unique <- Multi_Tp_sig %>% filter(!gene_id %in% (D59_w_multi %>% filter(Analysis_ID == a) %>% pull(gene_id)))
  multi_unique <- multi_unique[c(3,7)]
  multi_unique$Direction <- NA
  multi_unique$Analysis_ID <- a
  multi_unique$O59_Direction <- NA
  multi_unique$Comparison <- "Unique_to_mulit_tp"
  D59_w_multi <- rbind(D59_w_multi, multi_unique)
}

#Get designated trajectory specific counts for each analysis and the comparison calls
D59_w_multi$Trajectory <- gsub("-.*","",D59_w_multi$Official_ID)
D59_w_multi_count <- D59_w_multi %>% dplyr::count(Trajectory, Analysis_ID, Comparison)
D59_w_multi_count$Comparison <- factor(D59_w_multi_count$Comparison,levels = c("Unique_to_mulit_tp",
                                                                               "Unique_Old59",
                                                                               "Unique_OldAlt",
                                                                               "Direct_Call_Differ",
                                                                               "Direct_Call_Match"))

#Get label positions for graph
D59_w_multi_count <- D59_w_multi_count %>% arrange(Trajectory, Analysis_ID,desc(Comparison))

D59_w_multi_count$Percent <- NA  
D59_w_multi_count$Percent[D59_w_multi_count$Trajectory == "Complex"] <- D59_w_multi_count[D59_w_multi_count$Trajectory == "Complex", 4]/1010
D59_w_multi_count$Percent[D59_w_multi_count$Trajectory == "LinearDown"] <- D59_w_multi_count[D59_w_multi_count$Trajectory == "LinearDown", 4]/2591
D59_w_multi_count$Percent[D59_w_multi_count$Trajectory == "LinearUp"] <- D59_w_multi_count[D59_w_multi_count$Trajectory == "LinearUp", 4]/2536

#Get label positions
D59_w_multi_count$Percent_pos <- 0
D59_percent <- NULL
for (a in unique(D59_w_multi_count$Analysis_ID)) {
  for (t in unique(D59_w_multi_count$Trajectory)) {
    categories <- D59_w_multi_count %>% filter(Trajectory == t & Analysis_ID == a)
    sum <- 0
    for (i in 1:nrow(categories)) {
      if(i == 1){
        categories[i,6] <- categories[i,5]/2
        sum <- sum+categories[i,5]
      }else{
        categories[i,6] <- categories[i,5]/2+sum
        sum <- sum+categories[i,5]
      }
    }
    D59_percent <- rbind(D59_percent, categories)
  }
}

#Determine if the labels will need some nudging
D59_percent$nudge <- 0
D59_percent$nudge_y <- 0

for (i in 1:(nrow(D59_percent)-1)) {
  if(abs(D59_percent[i,6]-D59_percent[i+1,6]) < 0.025){
    D59_percent[i,7] <- 0.2
    D59_percent[i+1,7] <- -0.2
  }
  if(abs(D59_percent[i,6]-D59_percent[i+1,6]) < 0.015){
    D59_percent[i,8] <- -0.02
    D59_percent[i+1,8] <- 0.02
  }
}

#Create facet labels for each of the alternative analyses
facet_labs <- c("Old = Day 10 and 14",
                "Old = Day 14 and 17",
                "Old = Day 17 and 23",
                "Old = Day 23 and 27",
                "Old = Day 27 and 31",
                "Old = Day 31 and 36", 
                "Old = Day 36 and 38",
                "Old = Day 38 and 42",
                "Old = Day 42 and 48",
                "Old = Day 48 and 50",
                "Old = Day 50 and 55")
names(facet_labs) <- unique(D59_percent$Analysis_ID)

D59_compare_plot <- 
  D59_percent%>% ggplot()+
  geom_bar(aes(x=Trajectory, y=Percent, fill=Comparison), stat = "identity")+
  facet_wrap(~Analysis_ID, labeller = labeller(Analysis_ID=facet_labs))+
  scale_fill_manual(values = c("gray88",brewer.pal(8,"Paired")[c(1,3,5,7)]),
                    labels = c("Not identified in either Old Analysis",
                               "Old = Day59 Analysis Unique",
                               "Old = Alternative Analysis Unique",
                               "Shared Genes, Calls Differ",
                               "Shared Genes, Calls Match"))+
  geom_label(aes(label = paste(round(Percent,3)*100,"%",sep = ""), x=Trajectory, y=Percent_pos),
             label.size = 0, fill="white", alpha=0,size=2, nudge_x = D59_percent$nudge,nudge_y = D59_percent$nudge_y)+
  scale_y_continuous(labels = c("0%", "25%", "50%", "75%", "100%"))+
  xlab("Trajectory Designation")+
  #ggtitle("Alternative Old vs Day 59 Old Analyses")+
  theme_bw()+
  theme(legend.position = c(0.88,0.15),
        legend.text = element_text(size=5.25),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        axis.text.x = element_text(size=6))

ggsave("Plots/SupFig9_Comparison_w_D59.png", D59_compare_plot, height = 6, width = 6, units = "in")  
