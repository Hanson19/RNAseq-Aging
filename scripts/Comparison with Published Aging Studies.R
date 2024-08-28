#Introduction####
#Title: Comparison with Young v Old published Studies

#Purpose: Compare our identified genes with previous aging expression studies
#and see if we find the same genes and how do the directions called in 
#these studies compare with our linear expression trajectories

#Setup: Each main section begins with what packages are required and
#expected directory structure. Code through identifying genes shared with
#other papers, identifying our unique genes, and how to make supplemental 
#table K and Supplemental Figure 5. 

#Created:6/19/2024, KMH
#Last Edited: 8/28/2024, KMH

#Shared Genes####
#Identify shared genes between our genes and other Aging expression papers
#For each dataset we are going associate the valid IDs with the paper's genes
#and also call if the gene increased or decreased in expression

library(tidyverse)

#Code expects a certain directory structure with this specific files:
# Master Directory/
#   deseq_results/
#     All Analyses/
#       gene_color_id_ref.csv
#   Aging Analysis Results/
#     Comparison with Other Studies/
#       Sig Genes FB Gene ID Valid.txt
#       Bajgiran et al/
#         Bajgiran_FB_GeneIDs_Day1_v_Day60_filtered.txt
#         Bajgoram_DEG_fold_change_text.csv
#       Bordet et al/
#         Bordet_FB_GeneIDs_filtered.txt
#         Bordet_Supplementary Table S1.csv
#       Carnes et al/
#         Carnes_FB_GeneIDs_1week_v_5week_filtered.txt
#         Carnes_RNAseq_Male_Anovas_txt.csv
#       Giradot et al/
#         Giradot_FB_Genes_Day3_v_Day40_filtered.txt
#         Giradot_head_probe.csv
#       Highfill et al/
#         Highfill_FB_GeneIDs_Young_v_Old_head_filtered.txt
#         Highfill paper Addfile10 TableS6 Head.txt
#         Highfill_FB_GeneIDs_Young_v_Old_body_filtered.txt
#         Highfill paper Addfile9 TableS5 Body.txt
#       Lai et al/
#         Lai_S3_SigGenes_rmdoubles_txt.csv
#         Lai_FB_GeneIDs_old_v_young_filtered.txt
#       Landis et al/
#         Landis_SelectedGenes_from_Table_Txt.csv
#         Landis_FB_GeneIDs_old_v_young_filtered.txt
#       Zane et al/
#         Zane_FB_GeneIDs_Day20_v_Day40NS_filtered.txt
#         Zane_DEG_20NS_v_40NS_text.csv

shared_gene_numbers <- NULL

##Add Validated IDs to Identified Genes####
Valid_Ids <- read.table("Aging Analysis Results/Comparison with Other Studies/Sig Genes FB Gene ID Valid.txt", header = TRUE)
Gene_Ids <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)
Gene_Ids <- Gene_Ids[-c(1)]
colnames(Valid_Ids)[1] <- "gene_id"
Gene_Ids_Valid <- merge(Valid_Ids, Gene_Ids, by=c("gene_id"))
write.csv(Gene_Ids_Valid,"Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv")

#Each comparison's code is slightly different due to differences in data frame setup
#General steps are to add the validated gene IDs to the paper's df
#Determine if the paper called the gene's expression up or down
#Pull out the genes that are shared between our data and their using the validated IDs
#Record how many genes are shared and how many genes from the published paper had valid IDs. 

##Zane et al. 2023####
#Females, Full Body, Day20 vs Day40 Nonsmurf flies
Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]

#Add Valid Genes to Zane data
Zane_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Zane et al/Zane_FB_GeneIDs_Day20_v_Day40NS_filtered.txt", header = TRUE)
Zane_gene_stats <- read.csv("Aging Analysis Results/Comparison with Other Studies/Zane et al/Zane_DEG_20NS_v_40NS_text.csv", header = TRUE)
colnames(Zane_Valid_IDs)[1] <- "gene_id"
colnames(Zane_gene_stats)[1] <- "gene_id"
Zane_Valid_stats <- inner_join(Zane_Valid_IDs, Zane_gene_stats, by=c("gene_id"))
#515 genes with valid ids who have significant padj value

#Call if the gene increased or decreased in expression
#Positive fold change means it increased, negative means it decreased
Zane_Valid_stats$Direction <- NA
Zane_Valid_stats$Direction[Zane_Valid_stats$log2FoldChange > 0] <- "Up"
Zane_Valid_stats$Direction[Zane_Valid_stats$log2FoldChange < 0] <- "Down"

#Pull out genes shared with our data
Zane_shared <- inner_join(Zane_Valid_stats[c(2,9)],Sig_Genes,  by=c("validated_id"))
#366 genes are shared
write.csv(Zane_shared, "Aging Analysis Results/Comparison with Other Studies/Zane et al/Zane_Shared_Genes.csv")

#Add numbers about shared genes to shared_gene_numbers
shared_gene_numbers <- as.data.frame(rbind(shared_gene_numbers, cbind("Zane","2023",37822253,nrow(Zane_Valid_IDs), nrow(Zane_shared),
                                                                      nrow(Zane_shared)/nrow(Zane_Valid_IDs), "F",
                                                                      "Full_Body","Day20vDay40_Nonsmurf_Flies")))
##Bordet et al. 2021####
#Females, Full Body, Day2 v Day40
Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]

Bordet_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Bordet et al/Bordet_FB_GeneIDs_filtered.txt", header = TRUE)
Bordet_gene_stats <- read.csv("Aging Analysis Results/Comparison with Other Studies/Bordet et al/Bordet_Supplementary Table S1.csv", header = TRUE)
colnames(Bordet_Valid_IDs)[1] <- c("EntrezID")

#Match genes up with valid IDs
Bordet_valid_stats <- inner_join(Bordet_Valid_IDs, Bordet_gene_stats, by=c("EntrezID"))
#When look at the count things don't match Bordet gene stats is 1225 and Bordet Valid genes
#is 1183
db_genes <- Bordet_valid_stats %>% dplyr::count(EntrezID) %>% filter(n > 1) #34 genes whose ID appears more than once
Bordet_valid_stats %>% filter(EntrezID %in% db_genes$EntrezID)
#74 rows. Heres first 10:
#    EntrezID validated_id current_symbol      Name Young1 Young2 Young3  Old1  Old2  Old3
# 1     32780  FBgn0003380             Sh        Sh   5.05   6.90   6.26  4.25  4.61  4.47
# 2     32780  FBgn0003380             Sh        Sh   5.05   6.90   6.26  4.25  4.61  4.47
# 3     32780  FBgn0003380             Sh        Sh   5.05   6.90   6.26  4.25  4.61  4.47
# 4     33824  FBgn0001128          Gpdh1      Gpdh  12.31  12.54  12.29 10.81 11.18 11.11
# 5     33824  FBgn0001128          Gpdh1      Gpdh  12.31  12.54  12.29 10.81 11.18 11.11
# 6     34409  FBgn0032234            gny       gny  10.09  10.82  10.69 12.01 12.10 12.17
# 7     34409  FBgn0032234            gny       gny  10.09  10.82  10.69 12.01 12.10 12.17
# 8     34485  FBgn0032297        CG17124   CG17124   9.62   9.83   9.67  7.96  7.98  8.72
# 9     34485  FBgn0032297        CG17124   CG17124   9.62   9.83   9.67  7.96  7.98  8.72
# 10    34981  FBgn0020414          Idgf3     Idgf3   7.46   7.89   8.16 10.47 10.92 10.85

#When we look we see the counts are same so we are going to keep genes, but make each gene have one row
Bordet_valid_stats <- unique(Bordet_valid_stats)
#1183 to 1144
#Still have issue because we only have 1143 valid IDs, but 1144 in valid_stats
Bordet_valid_stats %>% dplyr::count(validated_id) %>% arrange(desc(n))
#     validated_id n
# 1    FBgn0263006 2
Bordet_valid_stats %>% filter(validated_id == "FBgn0263006")
#   EntrezID validated_id current_symbol    Name Young1 Young2 Young3  Old1  Old2  Old3
# 1    49297  FBgn0263006          SERCA Ca-P60A  11.46  11.52  11.52 10.45 10.06 10.42
# 2    49297  FBgn0263006          SERCA   SERCA  14.36  14.40  14.31 13.11 10.39 10.79

#We see that this gene id has two different symbols and the counts do not match each other.
#Remove this gene because we do not know which are the correct counts
Bordet_valid_stats <- Bordet_valid_stats %>% filter(validated_id != "FBgn0263006")
#1142 rows

#Determine if expression increased or decreased between old and young
Bordet_valid_stats$Change_Expression <- NA
for (i in 1:nrow(Bordet_valid_stats)) {
  Bordet_valid_stats[i,11] <- mean(Bordet_valid_stats[i,10], Bordet_valid_stats[i,9], Bordet_valid_stats[i,8])-
    mean(Bordet_valid_stats[i,5],Bordet_valid_stats[i,6], Bordet_valid_stats[i,7])
} 
Bordet_valid_stats$Direction[Bordet_valid_stats$Change_Expression > 0] <- c("Up")
Bordet_valid_stats$Direction[Bordet_valid_stats$Change_Expression < 0] <- c("Down")

#Pull out genes with shared with out data
Bordet_shared <- inner_join(Bordet_valid_stats[c(2,12)], Sig_Genes, by=c("validated_id"))
#690 genes shared between analyses
write.csv(Bordet_shared, "Aging Analysis Results/Comparison with Other Studies/Bordet et al/Bordet_Shared_Genes.csv")

#Add shared gene information to table
shared_gene_numbers <- rbind(shared_gene_numbers, cbind("Bordet","2021",34946931,nrow(Bordet_valid_stats), nrow(Bordet_shared),
                                                        nrow(Bordet_shared)/nrow(Bordet_valid_stats), "F",
                                                        "Full_Body", "Day2vDay45"))

##Bajgiran et al. 2021####
#M and F, Fullbody, D1vD60

Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]

#Add Validated IDs to gene stats
Bajgiran_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Bajgiran et al/Bajgiran_FB_GeneIDs_Day1_v_Day60_filtered.txt", header = TRUE)
Bajgiran_gene_stats <- read.csv("Aging Analysis Results/Comparison with Other Studies/Bajgiran et al/Bajgiran_DEG_fold_change_text.csv", header = TRUE)
colnames(Bajgiran_gene_stats)[1] <- "gene_id"
colnames(Bajgiran_Valid_IDs)[1] <- "gene_id"
Bajgiran_valid_stats <- inner_join(Bajgiran_Valid_IDs, Bajgiran_gene_stats, by=c("gene_id"))
#10,188 valid IDs

#Determine if expression increased or decreased using LogFC
Bajgiran_valid_stats$Direction[Bajgiran_valid_stats$LogFC > 0] <- "Up"
Bajgiran_valid_stats$Direction[Bajgiran_valid_stats$LogFC < 0] <- "Down"

#Pull out shared genes
Bajgiran_Shared <- inner_join(Bajgiran_valid_stats[c(2,13)], Sig_Genes, by=c("validated_id"))
#Share 4,814 genes
write.csv(Bajgiran_Shared, "Aging Analysis Results/Comparison with Other Studies/Bajgiran et al/Bajgiran_Shared_Genes.csv")

shared_gene_numbers <- rbind(shared_gene_numbers, cbind("Bajgiran", "2021", 34632013, nrow(Bajgiran_valid_stats), nrow(Bajgiran_Shared),
                                                        nrow(Bajgiran_Shared)/nrow(Bajgiran_valid_stats), "M&F", "Full_Body", "Day1vDay60"))

##Highfill et al. 2016####
#Has both data from heads and fully body
Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]
#Head
#Add validated IDs to gene stats
Highfill_head_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Highfill et al/Highfill_FB_GeneIDs_Young_v_Old_head_filtered.txt",  header = TRUE)
Highfill_head_gene_stats <- read.table("Aging Analysis Results/Comparison with Other Studies/Highfill et al/Highfill paper Addfile10 TableS6 Head.txt", header = TRUE)
colnames(Highfill_head_gene_stats)[2] <- "gene_id"
colnames(Highfill_head_Valid_IDs)[1] <- "gene_id"
Highfill_valid_stats <- inner_join(Highfill_head_Valid_IDs, Highfill_head_gene_stats, by=c("gene_id"))
#1,932 valid genes

#Call Direction expression change
Highfill_valid_stats$Direction[Highfill_valid_stats$FoldChange.Log2 > 0] <- c("Up")
Highfill_valid_stats$Direction[Highfill_valid_stats$FoldChange.Log2 < 0] <- c("Down")

#Pull out shared Genes
Highfill_head_shared <- inner_join(Highfill_valid_stats[c(2,14)], Sig_Genes, by=c("validated_id"))
#Shared 1,696 genes
write.csv(Highfill_head_shared, "Aging Analysis Results/Comparison with Other Studies/Highfill et al/Highfill_Shared_Genes_Head.csv")

shared_gene_numbers <- rbind(shared_gene_numbers, cbind("Highfill", "2016", 27485207, nrow(Highfill_valid_stats), nrow(Highfill_head_shared),
                                                        nrow(Highfill_head_shared)/nrow(Highfill_valid_stats), "F", "Head", "Day1:3vs50%Pop"))
#Body
#Add validated IDs to gene stats
Highfill_body_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Highfill et al/Highfill_FB_GeneIDs_Young_v_Old_body_filtered.txt", header = TRUE)
Highfill_body_gene_stats <- read.table("Aging Analysis Results/Comparison with Other Studies/Highfill et al/Highfill paper Addfile9 TableS5 Body.txt", header = TRUE)
colnames(Highfill_body_gene_stats)[2] <- "gene_id"
colnames(Highfill_body_Valid_IDs)[1] <- "gene_id"
Highfill_body_valid_stats <- inner_join(Highfill_body_Valid_IDs, Highfill_body_gene_stats, by=c("gene_id"))

#Have genes appearing multiple times. This is because they did several different comparisons
db_body <- Highfill_body_valid_stats %>% count(gene_id) %>% filter(n>1)
Highfill_body_valid_stats %>% filter(gene_id %in% db_body$gene_id)

#Call Direction expression change
Highfill_body_valid_stats$Direction[Highfill_body_valid_stats$FoldChange.Log2 > 0] <- c("Up")
Highfill_body_valid_stats$Direction[Highfill_body_valid_stats$FoldChange.Log2 < 0] <- c("Down")

Shared_Highfill_body <- inner_join(Highfill_body_valid_stats[c(2,14)], Sig_Genes, by=c("validated_id"))
#209 genes

#Solve duplicate gene
#If a gene's direction matches in the multiple cases keep gene
unique_shared_body <- unique(Shared_Highfill_body)
#158 genes

#Pull out genes they still have duplicates. These are ones where direction dind't
#match
df_body_2 <- unique_shared_body %>% dplyr::count(validated_id) %>% filter(n > 1)
#7 genes
#Remove 7 genes from our final set of genes
unique_shared_body <- unique_shared_body %>% filter(!validated_id %in% df_body_2$validated_id)
#144 genes
#9 genes whose direction is NA
gene_dir_na <- unique_shared_body %>% filter(is.na(Direction))
#These are all genes where one of the FPKM are 0
nas_og_Data <- Highfill_body_gene_stats %>% filter(gene_id %in% gene_dir_na$gene_id)
#Determine direction using FPKM
nas_og_Data$Direction[nas_og_Data$FPKM1 < nas_og_Data$FPKM2] <- "Up"
nas_og_Data$Direction[nas_og_Data$FPKM1 > nas_og_Data$FPKM2] <- "Down"
#Pull out gene id and the called direction 
nas_og_directed <- unique(nas_og_Data[,c(2,12)])
#Replace the direction nas when the newly determined directions
for (r in 1:nrow(unique_shared_body)) {
  if (is.na(unique_shared_body[r,2])){
    unique_shared_body[r,2] <- nas_og_directed %>% filter(gene_id == unique_shared_body[r,1]) %>% pull(Direction)
  }
}
unique_shared_body %>% filter(gene_id %in% nas_og_directed$gene_id)

write.csv(unique_shared_body, "Aging Analysis Results/Comparison with Other Studies/Highfill et al/Highfill_Shared_Genes_Body.csv")

#Subtracting 7 because we had to remove 7 genes from thhe pool
shared_gene_numbers <- rbind(shared_gene_numbers, cbind("Highfill", "2016", 27485207, nrow(Highfill_body_Valid_IDs)-7, nrow(unique_shared_body),
                                                        nrow(unique_shared_body)/(nrow(Highfill_body_Valid_IDs)-7), "F", "Body", "Day1:3vs50%Pop"))


#
##Carnes et al. 2015####
#M Full Body 1week vs 5 weeks

Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]

#Add validated IDs to gene stats
Carnes_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Carnes et al/Carnes_FB_GeneIDs_1week_v_5week_filtered.txt", header = TRUE)
Carnes_gene_stats <- read.csv("Aging Analysis Results/Comparison with Other Studies/Carnes et al/Carnes_RNAseq_Male_Anovas_txt.csv", header = TRUE)
colnames(Carnes_Valid_IDs)[1] <- "gene_id"
colnames(Carnes_gene_stats)[1] <- "gene_id"

Carnes_valid_stats <- inner_join(Carnes_Valid_IDs, Carnes_gene_stats, by=c("gene_id"))
#9,161 validated genes

#Determine if expression increased or decreased 
Carnes_valid_stats$Change <- Carnes_valid_stats$Week.5.Mean-Carnes_valid_stats$Week.1.Mean
Carnes_valid_stats$Direction[Carnes_valid_stats$Change > 0] <- c("Up")
Carnes_valid_stats$Direction[Carnes_valid_stats$Change < 0] <- c("Down")

Carnes_Shared <- inner_join(Carnes_valid_stats[c(2,26)], Sig_Genes, by=c("validated_id"))
#Share 4,308 genes
write.csv(Carnes_Shared, "Aging Analysis Results/Comparison with Other Studies/Carnes et al/Carnes_Shared_Genes.csv")

shared_gene_numbers <- rbind(shared_gene_numbers, cbind("Carnes", "2015", 26378456, nrow(Carnes_valid_stats), nrow(Carnes_Shared),
                                                        nrow(Carnes_Shared)/nrow(Carnes_valid_stats), "M", "Full_Body", "1weekv5weeks"))
##Lai et al. 2007####
Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]

Lai_gene_stats <- read.csv("Aging Analysis Results/Comparison with Other Studies/Lai et al/Lai_S3_SigGenes_rmdoubles_txt.csv", header = TRUE)
Lai_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Lai et al/Lai_FB_GeneIDs_old_v_young_filtered.txt", header = TRUE)
colnames(Lai_Valid_IDs)[1] <- "gene_id"
colnames(Lai_gene_stats)[7] <- "gene_id"

#Do have duplicate genes. Removing those without an id
Lai_gene_stats %>% dplyr::count(gene_id) %>% filter(n > 1)
Lai_gene_stats <- Lai_gene_stats %>% filter(gene_id != "---") 

#Add validated ids 
Lai_valid_stats <- inner_join(Lai_Valid_IDs, Lai_gene_stats, by=c("gene_id"))

#Call direction change
Lai_valid_stats$Direction[Lai_valid_stats$Fold.Change..Old.Young. > 1] <- c("Up")
Lai_valid_stats$Direction[Lai_valid_stats$Fold.Change..Old.Young. < 1] <- c("Down")
#Correct ones where direction is NA
Lai_valid_stats$Direction[is.na(Lai_valid_stats$Direction) & Lai_valid_stats$Old > Lai_valid_stats$Young] <- "Up"
Lai_valid_stats$Direction[is.na(Lai_valid_stats$Direction) & Lai_valid_stats$Old < Lai_valid_stats$Young] <- "Down"

#If the expression direction matches among duplicates we will keep gene
Lai_valid_stats<- unique(Lai_valid_stats[c(1:3,12)])
db_genes <- Lai_valid_stats %>% dplyr::count(gene_id) %>% filter(n>1)
#       gene_id n
# 1 FBgn0052486 2
Lai_valid_stats <- Lai_valid_stats %>% filter(gene_id != db_genes[1,1])
#2,113 validated IDs

Lai_shared <- inner_join(Lai_valid_stats[c(2,4)], Sig_Genes, by=c("validated_id"))
#share 1327 genes

write.csv(Lai_shared, "Aging Analysis Results/Comparison with Other Studies/Lai et al/Lai_Shared_Genes.csv")

shared_gene_numbers <- rbind(shared_gene_numbers, cbind("Lai", "2007", 17196240, nrow(Lai_valid_stats), nrow(Lai_shared),
                                                        nrow(Lai_shared)/nrow(Lai_valid_stats), "M&F", "Full_Body", "Day7v10%Pop"))
##Giradot et al. 2006####
Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]

Giradot_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Giradot et al/Giradot_FB_Genes_Day3_v_Day40_filtered.txt", header = TRUE)
Giradot_head_probe <- read.csv("Aging Analysis Results/Comparison with Other Studies/Giradot et al/Giradot_head_probe.csv", header = TRUE)
#Have modified from the original table so only showing head comparison and have ordered from smallest to greatest.

#Pull out probes significant for the head and mark as down or up
Giradot_Down <- Giradot_head_probe[1:1081,] #rows chosen based on coloring scheme in excel document.
Giradot_Up <- Giradot_head_probe[3269:4503,]
Giradot_Down$Direction <- "Down"
Giradot_Up$Direction <- "Up"
Giradot_Direction <- rbind(Giradot_Up, Giradot_Down)

#Match with validated IDs
colnames(Giradot_Valid_IDs)[1] <- "FlyBase_ID"
Giradot_Valid <- inner_join(Giradot_Direction, Giradot_Valid_IDs, by=c("FlyBase_ID"))
#2260 genes
#Have more rows than valid IDs. Check if genes with multiple probe sets have same direction
unique(Giradot_Valid[c(2,12,13)]) %>% dplyr::count(FlyBase_ID) %>% filter(n > 1)
# [1] FlyBase_ID n         
# <0 rows> (or 0-length row.names)
#So each row is agreeing in the direction so am keeping the genes

#Pull out shared genes
Giradot_shared <- inner_join(unique(Giradot_Valid[c(13,12)]), Sig_Genes, by=c("validated_id"))
#1489 shared 2222 Giradot had validated
write.csv(Giradot_shared, "Aging Analysis Results/Comparison with Other Studies/Giradot et al/Giradot_Shared_Genes.csv")

shared_gene_numbers <- as.data.frame(rbind(shared_gene_numbers, cbind("Giradot","2006",16584578,nrow(Giradot_Valid_IDs), nrow(Giradot_shared),
                                                                      nrow(Giradot_shared)/nrow(Giradot_Valid_IDs), "M",
                                                                      "Head","Day3vDay40")))
#
##Landis et al. 2004####
#M Full Body D10 vs D61  
Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]

#Add Validated Ids to stats
Landis_gene_stats <- read.csv("Aging Analysis Results/Comparison with Other Studies/Landis et al/Landis_SelectedGenes_from_Table_Txt.csv", header = TRUE)
Landis_gene_stats <- Landis_gene_stats[,c(1:11)]
colnames(Landis_gene_stats)[10] <- c("gene_id")
Landis_Valid_IDs <- read.table("Aging Analysis Results/Comparison with Other Studies/Landis et al/Landis_FB_GeneIDs_old_v_young_filtered.txt", header = TRUE)
colnames(Landis_Valid_IDs)[1] <- "gene_id"

Landis_valid_stats <- inner_join(Landis_Valid_IDs, Landis_gene_stats, by=c("gene_id"))
Landis_valid_stats %>% dplyr::count(gene_id) %>% filter(n > 1)
Landis_valid_stats %>% filter(gene_id == "FBgn0003357")

Landis_valid_stats$Direction[Landis_valid_stats$Old.Y > 0] <- c("Up")
Landis_valid_stats$Direction[Landis_valid_stats$Old.Y < 0] <- c("Down")
#Direction match up for duplicates so will keep
Landis_valid_stats[c(1,14)] %>% dplyr::count(gene_id, Direction) %>% filter(n > 1)
Landis_valid_stats <- unique(Landis_valid_stats[c(2,14)])
#848 validated genes

Landis_shared <- inner_join(Landis_valid_stats, Sig_Genes, by=c("validated_id"))
#402 shared genes

write.csv(Landis_shared, "Aging Analysis Results/Comparison with Other Studies/Landis et al/Landis_Shared_Genes.csv")

shared_gene_numbers <- rbind(shared_gene_numbers, cbind("Landis", "2004", 15136717, nrow(Landis_valid_stats), nrow(Landis_shared),
                                                        nrow(Landis_shared)/nrow(Landis_valid_stats), "M", "Full_Body", "Day10vDay61"))

##Sup Table L Partial####
#should have shared_gene_numbers table that has line for each paper
colnames(shared_gene_numbers) <- c("First_Author", "Year","PMID","Num_Valid_Genes", 
                                   "Num_Shared_Genes", "Percent_Shared", "Sex", 
                                   "Tissue", "Comparison")
write.csv(shared_gene_numbers,"Aging Analysis Results/Comparison with Other Studies/Shared_Genes_w_Papers.csv")

#The additional information of where the data was found was added manually
#after table was created

#Percent Genes Shared####
#Get number of our genes that are identified in the studies
library(tidyverse)
library(RColorBrewer)

#Master Directory/
# Aging Analysis Results/
#   Comparison with Other Studies/
#     First author et al/
#       First author_shared_Genes.csv

paper_compare_directory <- paste("Aging Analysis Results/Comparison with Other Studies")

Sig_Genes <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)
Sig_Genes <- Sig_Genes[-c(1)]

#Create a list to hold the found genes, and one to hold not found genes
Genes_Found <- NULL
Genes_Not_Found <- Sig_Genes

#For every directory in our comparison with other studies directory
for (DIR in list.dirs(paper_compare_directory)[-c(1)]) {
  print(DIR)
  #Make a list of dataframes that have "_shared_"in that direcotry.
  shared_lists <- list.files(path = paste(DIR,"/", sep = ""), pattern = "_shared_", ignore.case = TRUE)
  if (length(shared_lists) != 0){
    #For each data frame in that list
    for (t in 1:length(shared_lists)) {
      #Read in the dataframes and pull out validated_id and Official cluster ID
      shared_genes_df <- read.csv(paste(DIR, "/", shared_lists[t], sep = ""), header = TRUE)
      shared_genes <- shared_genes_df[c("validated_id", "Official_ID")]
      #In our Genes Not Found, filter out all genes that are in the shared_genes datframe
      Genes_Not_Found <- Genes_Not_Found %>% filter(!validated_id %in% shared_genes$validated_id)
      #Add shared_genes along with the data frames name to Genes_Found
      shared_genes <- as.data.frame(cbind(shared_genes, shared_lists[t]))
      colnames(shared_genes) <- c("validated_id","Official_ID","file")
      Genes_Found <- rbind(Genes_Found, shared_genes)
    }
  }
}

dim(Genes_Found)
# [1] 15236     3
dim(Genes_Not_Found)
# [1] 287   9
#287 genes that are unique to our analysis

#Save Tables
write.csv(Genes_Found, "Aging Analysis Results/Comparison with Other Studies/Shared_Genes.csv")
write.csv(Genes_Not_Found, "Aging Analysis Results/Comparison with Other Studies/Unique_Genes.csv")                                               

###Sup Figure 11####
#Make bar graph for each of our designations that shows how many times their genes
#were identified in the 9 dataframes

library(tidyverse)
library(RColorBrewer)
Genes_Found <- read.csv("Aging Analysis Results/Comparison with Other Studies/Shared_Genes.csv", header = TRUE)
Genes_Not_Found <- read.csv("Aging Analysis Results/Comparison with Other Studies/Unique_Genes.csv", header = TRUE)
Genes_Found_Tally <- Genes_Found %>% dplyr::count(validated_id, Official_ID) %>% arrange(desc(n))

Not_Found_Tally <- as.data.frame(cbind(Genes_Not_Found$validated_id, Genes_Not_Found$Official_ID, 0))
colnames(Not_Found_Tally) <- c("validated_id","Official_ID", "n")

#Combine dataframes
Gene_Tally <- rbind(Genes_Found_Tally, Not_Found_Tally)
#Create a designation column based off the official IDs
Gene_Tally$Designation <- gsub("-.*","", Gene_Tally$Official_ID)
Gene_Tally$Count <- as.character(Gene_Tally$n)
#Get the gene identification counts and Total gene counts for trajectories
Designation_Count <- Gene_Tally %>% dplyr::count(Designation, Count)
Designation_TG <- Gene_Tally %>% dplyr::count(Designation)
colnames(Designation_TG)[2] <- "Total_Genes"
#Merge and then get the percent composition for each trajectory designation
Designation_Count <- merge(Designation_Count, Designation_TG, by=c("Designation"))  
Designation_Count$Percent <- Designation_Count$n/Designation_Count$Total_Genes*100

Identification_Trajectory <- 
  Designation_Count %>% ggplot(aes(x=Designation, y=Percent, fill=Count))+
  geom_bar(stat = "identity")+
  ggtitle("Gene Identification in Previous Papers")+
  scale_fill_manual(values =c("darkmagenta","gray98",brewer.pal(9, "Greys")[2:9]), name="Times\nFound")+
  scale_y_continuous(labels = c("0%", "25%", "50%", "75%", "100%"))+
  theme_bw()

ggsave("Plots/SupFig11_Identification_Trajectory.png", Identification_Trajectory, height = 6, width = 6, units = "in")

designation_aov <- aov(lm(Percent~Designation*Count, data = Designation_Count))
anova(designation_aov)
# Response: Percent
#                   Df Sum Sq Mean Sq F value Pr(>F)
# Designation        2    7.9    3.97     NaN    NaN
# Count              9 3251.9  361.32     NaN    NaN
# Designation:Count 16   38.1    2.38     NaN    NaN
# Residuals          0    0.0     NaN  

#Direction Counts####
#Supp Figure 12####
#For our shared genes, get the direction counts based on trajectory designations

library(tidyverse)
#Master Directory/
# Aging Analysis Results/
#   Comparison with Other Studies/
#     First author et al/
#       First author_shared_Genes.csv

direction_match <- function(DIR, NIX=NULL){
  #Make a dataframe that combines most of the information in each _shared_Genes.csv
  #dataframe 
  master_direction_list <- NULL
  for (d in list.dirs(DIR)[-c(1)]) {
    shared_gene_df <- list.files(path = paste(d, "/", sep = ""), pattern = "_shared_", ignore.case = TRUE)
    if (length(shared_gene_df) != 0){
      for (df in 1:length(shared_gene_df)) {
        if (!shared_gene_df[df] %in% NIX){
          shared_genes <- read.csv(paste(d, "/",shared_gene_df[df], sep = ""), header = TRUE,colClasses = c("character"))
          shared_gene_extract <- shared_genes[,c("validated_id","current_symbol","Direction","Official_ID","Cluster", "Dendo_clus", "col",
                                                 "New_color", "Analysis_ID")]
          shared_gene_extract$df <- shared_gene_df[df]
          master_direction_list <- rbind(master_direction_list, shared_gene_extract)
        }
      }
    }
  }
  #Give each gene their trajectory designation using first part of their cluster's ID
  master_direction_list$Designation <- gsub("-.*","", master_direction_list$Official_ID)
  
  #Count how many genes shared with paper, are in each trajectory designation and called either up or down
  master_direction_count <- master_direction_list %>% dplyr::count(Designation,Direction,df)
  master_direction_count <- master_direction_count %>% filter(!is.na(Direction))
  
  #Read in our summary shared gene table
  shared_df <- read.csv("Aging Analysis Results/Comparison with Other Studies/Shared_Genes_w_Papers.csv", header = TRUE)
  
  #Create a first author column in master direction count and then pull out information regarding
  #sex, tissue , and what was the comparison from shared_df
  master_direction_count$First_Author <- gsub("_.*","", master_direction_count$df)
  for (a in unique(master_direction_count$First_Author)) {
    print(a)
    if(a != "Highfill"){
      master_direction_count$Sex[master_direction_count$First_Author == a] <- shared_df %>% filter(First_Author == a) %>% pull(Sex)
      master_direction_count$Tissue[master_direction_count$First_Author == a] <- shared_df %>% filter(First_Author == a) %>% pull(Tissue)
      master_direction_count$Analysis[master_direction_count$First_Author == a] <- shared_df %>% filter(First_Author == a) %>% pull(Comparison)
      master_direction_count$y_post[master_direction_count$First_Author == a] <- max(master_direction_count %>% filter(First_Author == a) %>% pull(n))
    }else{ #High needs own section due to using 2 tissues
      master_direction_count$Sex[master_direction_count$df == "Highfill_Shared_Genes_Body.csv"] <- shared_df %>% filter(First_Author == a & Tissue == "Body") %>% pull(Sex)
      master_direction_count$Tissue[master_direction_count$df == "Highfill_Shared_Genes_Body.csv"] <- shared_df %>% filter(First_Author == a & Tissue == "Body") %>% pull(Tissue)
      master_direction_count$Analysis[master_direction_count$df == "Highfill_Shared_Genes_Body.csv"] <- shared_df %>% filter(First_Author == a & Tissue == "Body") %>% pull(Comparison)
      master_direction_count$y_post[master_direction_count$df == "Highfill_Shared_Genes_Body.csv"] <- max(master_direction_count %>% filter(df == "Highfill_Shared_Genes_Body.csv") %>% pull(n))
      
      master_direction_count$Sex[master_direction_count$df == "Highfill_Shared_Genes_Head.csv"] <- shared_df %>% filter(First_Author == a & Tissue == "Head") %>% pull(Sex)
      master_direction_count$Tissue[master_direction_count$df == "Highfill_Shared_Genes_Head.csv"] <- shared_df %>% filter(First_Author == a & Tissue == "Head") %>% pull(Tissue)
      master_direction_count$Analysis[master_direction_count$df == "Highfill_Shared_Genes_Head.csv"] <- shared_df %>% filter(First_Author == a & Tissue == "Head") %>% pull(Comparison)
      master_direction_count$y_post[master_direction_count$df == "Highfill_Shared_Genes_Head.csv"] <- max(master_direction_count %>% filter(df == "Highfill_Shared_Genes_Head.csv") %>% pull(n))
    }
  }
  #Create labels for each shared papers gene and clean up some names by removing _s
  master_direction_count$y_post <- master_direction_count$y_post-(.1*master_direction_count$y_post)
  master_direction_count$Analysis[master_direction_count$First_Author == "Zane"] <- gsub("_Non.*","", 
                                                                                         unique(master_direction_count %>% filter(First_Author == "Zane") %>% pull(Analysis)))
  master_direction_count$Tissue <- gsub("_"," ", master_direction_count$Tissue)
  master_direction_count$Title <- gsub(".csv", "", master_direction_count$df)
  master_direction_count$Title <- gsub("_", " ", master_direction_count$Title)
  
  master_direction_count$Labels <- paste(master_direction_count$Sex," ",master_direction_count$Tissue,"\n",master_direction_count$Analysis, sep = "")
  
  direction_plot <- master_direction_count %>% ggplot(aes(x=Designation, y=n, fill=Direction))+
    geom_bar(stat = "identity", position = "dodge")+
    facet_wrap(~Title, scales = "free")+
    xlab("Trajectory Designation")+
    ylab("Number of Genes")+
    geom_label(aes(label=Labels, y=y_post),x="Complex",fill="white")+
    #ggtitle("Old v Young Expression Change Designation")+
    scale_fill_manual(values = c("tomato", "skyblue"), name="Young v. Old Analysis Call")+
    theme_bw()+
    theme(legend.position = "bottom", text = element_text(size = 15))
  ggsave("Aging Analysis Results/Comparison with Other Studies/Trajectory_Direction_counts.png",direction_plot,width = 12, height = 12, units = "in")
  write.csv(master_direction_count,"Aging Analysis Results/Comparison with Other Studies/Trajectory_Direction_counts.csv")
}

direction_match(DIR = paste("Aging Analysis Results/Comparison with Other Studies"), NIX = NULL)

#All Identified Genes####
library(tidyverse)
#Master Directory/
# Aging Analysis Results/
#   Comparison with Other Studies/
#     First author et al/
#       First author_FB_GeneIDs_filtered.txt

#This function takes all the Validated IDs across 9 dataframes and list them
#together, along with counting how many time each gene is identified
comparison_study_genes <- function(DIR, NIX=NULL){
  master_gene_list <- NULL
  for (d in list.dirs(DIR)[-c(1)]) {
    if (!d %in% NIX){
      validated_df <- list.files(path = paste(d,"/",sep = ""), pattern = "_filtered", ignore.case = TRUE)
      if (length(validated_df) != 0){
        for (df in 1:length(validated_df)) {
          identified_genes <- read.table(paste(d, "/", validated_df[df],sep = ""), header = TRUE, colClasses = "character")
          valid_id <- identified_genes$validated_id
          for (id in 1:length(valid_id)) {
            if(!valid_id[id] %in% rownames(master_gene_list)){
              row_df <- as.data.frame(cbind(valid_id[id], 1))
              row_df <- column_to_rownames(row_df, "V1")
              master_gene_list <- rbind(master_gene_list, row_df)
            }else{
              master_gene_list$V2 <- as.numeric(master_gene_list$V2)
              master_gene_list[valid_id[id],1] <- master_gene_list[valid_id[id],1]+1
            }
          }
        }
      }
    }
  }

  master_gene_list <- rownames_to_column(master_gene_list, "gene_id")
  colnames(master_gene_list)[2] <- "Times_ID"
  print(master_gene_list)
}

gene_list <- comparison_study_genes(DIR = paste("Aging Analysis Results/Comparison with Other Studies"))
gene_list %>% dplyr::count(Times_ID)
#   Times_ID    n
# 1        1 3693
# 2        2 5069
# 3        3 2170
# 4        4  995
# 5        5  429
# 6        6  187
# 7        7   75
# 8        8   27
# 9        9    5
nrow(gene_list)
#12650 across 9 dataframes that had valid IDs

Sig_Genes_multitimepoint <- read.csv("Aging Analysis Results/Comparison with Other Studies/Sig_Genes_Valid_IDs.csv", header = TRUE)

identified_genes <- gene_list %>% filter(gene_id %in% Sig_Genes_multitimepoint$validated_id)
nrow(identified_genes)
#We identified 5855 of them
