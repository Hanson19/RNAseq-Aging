#Introduction####
#Title: Cluster Enrichment Analysis 

#Purpose: This script will walk through how we prepared the two lists needed
#to run the enrichment analysis using PANGEA 
#(https://www.flyrnai.org/tools/pangea/web/multiple_search/7227)
#and the three sets of code that seperated significant terms based on clusters,
#unique terms found in our three trajectory designations (LinearUp,LinearDown, 
#and Complex), and unique terms found after the first divide in dendrogram
#which is primarily based on if trajectories go up or down. 

#Setup:
#Each section begins with brief description of the goal, what R packages are 
#required and the directory structure including files the code expects. 
#Code does require 3 files made in "Identifying_and_Clustering.R":
#     day_NormCounts_deseq.csv
#     gene_color_id_ref.csv
#     cluster_color_id_ref.csv
#     Combo_Analysis_Clustered_Z_scores_k28.csv

#Created: 6/12/2024, KMH
#Last Edited: 11/06/2024, KMH

#PANGEA Preparation####
#We will be using PANGEA (https://www.flyrnai.org/tools/pangea/web/multiple_search/7227)
#multiple gene list search to identify terms enriched in our clusters. 
#We will also be doing our enrichment analysis with a background
#We need to get list of our significant genes with their cluster listed
#and we need to get list of genes for our background. 

#Packages
library(tidyverse)

#Code expects a certain directory structure with these specific files:
# Master Directory/
#   deseq_results/
#     day/
#       day_NormCounts_deseq.csv
#     All Analyses/
#       gene_color_id_ref.csv

dir.create("Aging Analysis Results")
dir.create("Aging Analysis Results/Enrichment")

#Identified Genes with Cluster list
gene_ids <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)
genes <- gene_ids$gene_id
clusters <- gene_ids$Official_ID
gene_cluster <- cbind(clusters, genes)
write.csv(gene_cluster,"Aging Analysis Results/Enrichment/pangea_id_and_cluster.csv", row.names = FALSE, col.names = FALSE)
#Manually go in and delete the first row which will saw genes, clusters
#Then open the file in Notepad or equivalent, copy file, and then paste into PANGEA

#Background genes
#Background will be all the genes which expression was detected in a sample
norm_counts <- read.csv("deseq_results/day/day_NormCounts_deseq.csv", header = TRUE)
norm_counts$sum <- NA
for (i in 1:nrow(norm_counts)) {
  norm_counts[i,50] <- sum(norm_counts[i,2:49])
}
exp_counts <- norm_counts %>% filter(sum != 0)
exp_genes <- exp_counts[1]
write.csv(exp_genes,"Aging Analysis Results/Enrichment/pangea_background_genes.csv", row.names=FALSE)
#Manually go in and delete the X that is the first row. 

#PANGEA Enrichment####
#When both the gene list with clusters and background gene list are made go to  
#https://www.flyrnai.org/tools/pangea/web/multiple_search/7227
#copy "pangea_id_and_cluster.csv" and paste into box 2 "Enter Genes and Groups" and paste
#gene ids in "pangea_background_genes" into box 3 "Enter customized Background". 

#We did 6 separate enrichment studies:
#SLIM2 GO BP
#SLIM2 GO CC
#SLIM2 GO MF
#DRSC GLAD Gene Group
#FlyBase Gene Group
#REACTOME pathway
#and downloaded the tables as csvs

#Going to pull out the enriched terms from these studies for each of the clusters.
#The code requires a specific directory setup.
#Each csv files is one of th files from PANGEA.

# Master Directory/
#   Aging Analysis Results/
#     Enrichment/
#       gene groups/
#         enrich_genegroups_DRSC_GLAD.csv 
#         enrich_genegroups_Flybase.csv
#       gene ontology/
#         enrich_SLIM2_BP.csv
#         enrich_SLIM2_CC.csv
#         enrich_SLIM2_MF.csv
#       pathways/
#         enrich_pathway_RECTOME.csv

library(tidyverse)
library(UpSetR)

pangea_enrichment_DeseqExp_Background <- function(GROUPING, SUBGROUP = NA, PVALVER = NULL, THRESHOLD = 0.05){
  #Read in the the appropriate pangea table and dataset
  allsig_genes <- read.csv("deseq_results/All Analyses/gene_color_id_ref.csv", header = TRUE)
  if (GROUPING == "Gene_Group" & SUBGROUP == "Flybase"){
    enrich_terms <- read.csv("Aging Analysis Results/Enrichment/gene groups/enrich_genegroups_Flybase.csv",header = TRUE)
  }else if (GROUPING == "Pathway" & SUBGROUP == "Reactome"){
    enrich_terms <- read.csv("Aging Analysis Results/Enrichment/pathways/enrich_pathway_REACTOME.csv", header = TRUE)
  }else if (GROUPING == "Ontology" & SUBGROUP == "SLIM2BP"){
    enrich_terms <- read.csv("Aging Analysis Results/Enrichment/gene ontology/enrich_SLIM2_BP.csv", header = TRUE)
  }else if (GROUPING == "Ontology" & SUBGROUP == "SLIM2MF"){
    enrich_terms <- read.csv("Aging Analysis Results/Enrichment/gene ontology/enrich_SLIM2_MF.csv", header = TRUE)
  }else if (GROUPING == "Ontology" & SUBGROUP == "SLIM2CC"){
    enrich_terms <- read.csv("Aging Analysis Results/Enrichment/gene ontology/enrich_SLIM2_CC.csv", header = TRUE)
  }else if (GROUPING == "Gene_Group" & SUBGROUP == "DRSC_GLAD"){
    enrich_terms <- read.csv("Aging Analysis Results/Enrichment/gene groups/enrich_genegroups_DRSC_GLAD.csv",header = TRUE)
  }
  #For every cluster
  for (c in str_replace(unique(allsig_genes$Official_ID), "-", ".")) {
    print(c)
    #Pull out columns in PANGEA data set that end with the cluster's ID
    clust_pangea <- enrich_terms[,str_ends(colnames(enrich_terms), c)]
    groups <- enrich_terms[,2:4]
    clust_pangea <- cbind(groups, clust_pangea)
    #We use are using Benjamini Hochberg false disocery rate. Designate that
    #column as the source of our p value
    if (PVALVER == "Benjamini_Hochberg_"){
      p_col_name <- paste("p.value..Benjamini...Hochberg..",c,sep = "")
    }else{
      p_col_name <- paste("p.value.",c,sep = "")
    }
    #pull out terms for that cluster whose p-value is less than 0.05
    sig_pangea <- clust_pangea[clust_pangea[p_col_name] <= THRESHOLD,]
    #This will also pull out NAs so removing those
    sig_pangea <- sig_pangea[!is.na(sig_pangea[p_col_name]),]
    #Reorinate pangea data sheet
    term_gene_long <- NULL
    speciesID_col_name <- paste("species.specific.ID.",c,sep = "")
    genes_group_wide <- NULL
    print("genes and grouping")
    #If there are significant terms, for very term
    if (nrow(sig_pangea) != 0){
      for (t in unique(sig_pangea$Gene.Set.ID)) {
        #Associate terms with genes
        #Pull out term specific information
        term_info <- unique(sig_pangea[sig_pangea["Gene.Set.ID"] == t, 1:3])
        #Get a list of the genes who match these term
        gene_vector <- str_split_1(sig_pangea %>% filter(Gene.Set.ID == t) %>% pull(speciesID_col_name),
                                   pattern = "; ")
        #Make list of genes into a dataframe and then create a column whose values are all 1
        gene_vector_wide <- as.data.frame(gene_vector)
        gene_vector_wide$Inclusion <- 1
        #Rename Inclusion column with the name of the gene set
        colnames(gene_vector_wide)[2] <- unique(sig_pangea %>% filter(Gene.Set.ID == t) %>% pull(Gene.Set.Name))
        #Combine this term speicifc list with data table with all terms
        if(is.null(genes_group_wide)){
          genes_group_wide <- gene_vector_wide
        }else{
          genes_group_wide <- full_join(genes_group_wide, gene_vector_wide, by=c("gene_vector"))
        }
        #Create a long data frame
        #For every gene associated with the term
        for (g in 1:length(gene_vector)) {
          #make a row with term info, gene name, and cluster
          gene_and_term_info <- cbind(term_info,gene_vector[g],c)
          #add row to term_ gene_long dataframe
          term_gene_long <- rbind(gene_and_term_info, term_gene_long)
        }
      }
      #There are genes which are not associated with terms and we want to designate them
      #pull out gene ids of genes that are not found in the gene groups wide dataframe
      not_found <- allsig_genes %>% filter(Official_ID == str_replace(c,"\\.","-")) %>% filter(!gene_id %in% genes_group_wide$gene_vector) %>% pull(gene_id)
      #Make this a data frame and give a column called not found with value of 1
      not_found <- as.data.frame(not_found)
      not_found$Not_Found_Category <- 1
      colnames(not_found)[1] <- c("gene_vector")
      #merge not found with gene_group wide
      genes_group_wide <- full_join(genes_group_wide, not_found, by=c("gene_vector"))
      #Anywhere in gene groups wide where there is a NA replace with a 0 
      genes_group_wide <- genes_group_wide %>% replace(is.na(.),0)
      #We can now use this dataframe for our upset graphs plotted below
      
      #Save our 3 data frame. Sig pangea terms, genes_group_wide, and term_gene_long
      print("All genes grouped now saving")
      #Alrighty so now we need to get the address which will differ depending
      if (GROUPING == "Gene_Group"){
        address <- paste("Aging Analysis Results/Enrichment/gene groups/",c,sep = "")
      }else if (GROUPING == "Ontology"){
        address <- paste("Aging Analysis Results/Enrichment/gene ontology/",c,sep = "")
      }else if (GROUPING == "Pathway"){
        address <- paste("Aging Analysis Results/Enrichment/pathways/",c,sep = "")
      }
      if (length(dir(address)) == 0){
        dir.create(address)
      }
      if (GROUPING == "Gene_Group"){
        write.csv(term_gene_long, paste(address,"/",PVALVER,"enrich_genegroups_",SUBGROUP,"_long_",c,".csv", sep = ""))
        write.csv(sig_pangea, paste(address,"/",PVALVER,"enrich_genegroups_",SUBGROUP,"_pangea_",c,".csv", sep = ""))
        write.csv(genes_group_wide, paste(address, "/",PVALVER,"enrich_genegroups_",SUBGROUP,"_wide_",c,".csv",sep = ""))
        pdf_title <- paste(address, "/",PVALVER,"enrich_genegroups_",SUBGROUP,"_upset_",c,".pdf", sep = "")
      }else if (GROUPING == "Ontology"){
        write.csv(term_gene_long, paste(address,"/", PVALVER,"enrich_ontology_",SUBGROUP,"_long_",c,".csv", sep = ""))
        write.csv(sig_pangea, paste(address, "/",PVALVER,"enrich_ontology_",SUBGROUP,"_pangea_",c,".csv",sep = ""))
        write.csv(genes_group_wide, paste(address, "/",PVALVER,"enrich_ontology_",SUBGROUP,"_wide_",c,".csv",sep = ""))
        pdf_title <- paste(address, "/",PVALVER,"enrich_ontology_",SUBGROUP,"_upset_",c,".pdf", sep = "")
      }else if (GROUPING == "Pathway"){
        write.csv(term_gene_long, paste(address, "/",PVALVER,"enrich_pathway_",SUBGROUP,"_long_",c,".csv",sep = ""))
        write.csv(sig_pangea, paste(address, "/",PVALVER,"enrich_pathway_",SUBGROUP,"_pangea_",c,".csv",sep = ""))
        write.csv(genes_group_wide, paste(address, "/",PVALVER,"enrich_pathway_",SUBGROUP,"_wide_",c,".csv",sep = ""))
        pdf_title <- paste(address, "/",PVALVER,"enrich_pathway_",SUBGROUP,"_upset_",c,".pdf", sep = "")
      }
      print("Upsetting")
      #Plot to see the overlap of significant terms with genes
      #make gene vector column into rownames of gene_group_wide
      genes_group_wide <- column_to_rownames(genes_group_wide, "gene_vector")
      #Create upset plot and save
      upsetdiag <- 
        UpSetR::upset(genes_group_wide, nsets = length(unique(sig_pangea$Gene.Set.ID))+1, number.angles = 30, point.size = 3.5, line.size = 2, 
                      mainbar.y.label = "Group Intersections", sets.x.label = "Genes Per Group", 
                      text.scale = c(1.75, 1.75, 1.5, 1.5, 1, 2))
      pdf(file = pdf_title, width = 15, height = 15)
      print(upsetdiag)
      graphics.off()
    }else{
      #If there are no significant terms, save a blank dataframe in cluster's folder
      print("No enriched groups")
      place_holder <- clust_pangea[1,]
      place_holder[1,] <- NA
      if (GROUPING == "Gene_Group"){
        address <- paste("Aging Analysis Results/Enrichment/gene groups/",c,sep = "")
        file_name <- paste("/",PVALVER,"enrich_genegroups_",SUBGROUP,"_pangea_",c,".csv", sep = "")
      }else if (GROUPING == "Ontology"){
        address <- paste("Aging Analysis Results/Enrichment/gene ontology/",c,sep = "")
        file_name <- paste("/",PVALVER,"enrich_ontology_", SUBGROUP,"_pangea_",c,".csv", sep = "")
      }else if (GROUPING == "Pathway"){
        address <- paste("Aging Analysis Results/Enrichment/pathways/",c,sep = "")
        file_name <- paste("/",PVALVER,"enrich_pathway_",SUBGROUP,"_pangea_",c,".csv",sep = "")
      }
      if (length(dir(address)) == 0){
        dir.create(address)
      }
      write.csv(place_holder, paste(address, file_name, sep = ""))
    }
  }
}

ben_hoch <- "Benjamini_Hochberg_"

#Grouping
gg <- "Gene_Group"
ont <- "Ontology"
path <- "Pathway"


#Subgroups
fb <- "Flybase"
slim2bp <- "SLIM2BP"
slim2mf <- "SLIM2MF"
slim2cc <- "SLIM2CC"
drsc <- "DRSC_GLAD"
rct <- "Reactome"

pangea_enrichment_DeseqExp_Background(GROUPING = ont, PVALVER = ben_hoch, SUBGROUP = slim2bp)
pangea_enrichment_DeseqExp_Background(GROUPING = ont, PVALVER = ben_hoch, SUBGROUP = slim2mf)
pangea_enrichment_DeseqExp_Background(GROUPING = ont, PVALVER = ben_hoch, SUBGROUP = slim2cc)

pangea_enrichment_DeseqExp_Background(GROUPING = gg, PVALVER = ben_hoch, SUBGROUP = fb)
pangea_enrichment_DeseqExp_Background(GROUPING = gg, PVALVER = ben_hoch, SUBGROUP = drsc)

pangea_enrichment_DeseqExp_Background(GROUPING = path, PVALVER = ben_hoch, SUBGROUP = rct)

#Trajectory Classification Divided####
#Identify unique terms found in our Complex, LinearUp and LinearDown Clusters
#Unlike PANGEA Enrichment section we can test all the enrichment classes at once

#The code requires a specific directory setup.
# Master Directory/
#   deseq_results/
#     All_Analyses/
#       cluster_color_id_ref.csv
#   Aging Analysis Results/
#     Enrichment/
#       gene groups/
#         Cluster.#/
#           Benjamini_Hochberg_enrich_genegroups_DRSC_GLAD_pangea_Cluster.#.csv
#           Benjamini_Hochberg_enrich_genegroups_Flybase_pangea_Cluster.#.csv
#       gene ontology/
#         Cluster.#/
#           Benjamini_Hochberg_enrich_ontology_SLIM2BP_pangea_Cluster.#.csv
#           Benjamini_Hochberg_enrich_ontology_SLIM2CC_pangea_Cluster.#.csv
#           Benjamini_Hochberg_enrich_ontology_SLIM2MF_pangea_Cluster.#.csv
#       pathways/
#         Cluster.#/
#           Benhamini_Hochberg_enrich_pathway_Reactome_pangea_Cluster.#.csv

library(tidyverse)

sorted_by_background <- function(GROUPING, SUBGROUPING){
  #Read in cluster id information
  clust_desg <- read.csv("deseq_results/All Analyses/cluster_color_id_ref.csv", header = TRUE)
  #Create new column just has the designation of Complex, LinearUp or Linear Down
  clust_desg$Designation_Direction <-gsub("-.*","",clust_desg$Official_ID)
  #Pull out the cluster's Ids that are in each of these trajectory designations
  linear_up <- clust_desg %>% filter(Designation_Direction == "LinearUp") %>% pull(Official_ID)
  linear_down <- clust_desg %>% filter(Designation_Direction == "LinearDown") %>% pull(Official_ID)
  complex <- clust_desg %>% filter(Designation_Direction == "Complex") %>% pull(Official_ID)
  #Create a master list to store the divided terms
  Divided_terms <- NULL
  #For each Enrichment term class
  for (g in GROUPING) {
    print(g)
    for (s in SUBGROUPING) {
      print(s)
      #Create list to store enriched_Terms for that specific class
      enrich_terms <- NULL
      #For each cluster
      for (c in str_replace(clust_desg$Official_ID, "-", ".")) { #Note at times need to go back and forth between Cluster-# and Cluster.#
        if(g == "Gene_Group"){ 
          cluster_term_add <- paste("Aging Analysis Results/Enrichment/gene groups/",c,"/Benjamini_Hochberg_enrich_genegroups_",
                                    s,"_pangea_",c,".csv",sep = "")
        }else if (g == "Ontology"){
          cluster_term_add <- paste("Aging Analysis Results/Enrichment/gene ontology/",c,"/Benjamini_Hochberg_enrich_ontology_",
                                    s,"_pangea_",c,".csv",sep = "")
        }else if (g == "Pathway"){
          cluster_term_add <- paste("Aging Analysis Results/Enrichment/pathways/",c,"/Benjamini_Hochberg_enrich_pathway_",
                                    s,"_pangea_",c,".csv",sep = "")
        }
        if (file.exists(cluster_term_add)){
          #Read in the Cluster's Significant PANGEA Terms
          cluster_term <- read.csv(cluster_term_add, header = TRUE)
          #Creat a row with the term information and name of the Cluster
          row <- cbind(cluster_term[,c(2:4)],str_replace(c,"\\.","-"))
          colnames(row)[4] <- c("Cluster")
          #Bind the rows to the class specific list
          enrich_terms <- rbind(enrich_terms, row)
        }
      }
      #If there are terms in enrich_Terms
      if(!is.null(enrich_terms)){
        #Divided terms based on which trajectory designation the cluster belongs to
        linup <- enrich_terms %>% filter(Cluster %in% linear_up)
        lindown <- enrich_terms %>% filter(Cluster %in% linear_down)
        comp <- enrich_terms %>% filter(Cluster %in% complex)
        
        #For each trajectory designation pull our the terms that are only unique
        #to that trajectory designation
        linup_unique <- linup %>% filter(!Gene.Set.ID %in% lindown$Gene.Set.ID &
                                           !Gene.Set.ID %in% comp$Gene.Set.ID)
        lindown_unique <- lindown %>% filter(!Gene.Set.ID %in% linup$Gene.Set.ID &
                                               !Gene.Set.ID %in% comp$Gene.Set.ID)
        complex_unique <- comp %>% filter(!Gene.Set.ID %in% lindown$Gene.Set.ID &
                                            !Gene.Set.ID %in% linup$Gene.Set.ID)
        #Create a list of the clusters that have the unique term.
        #Will do this for each trajectory so only comments on LinearUp's loop
        for (t in unique(linup_unique$Gene.Set.ID)) {
          term <- linup_unique %>% filter(Gene.Set.ID == t)
          for (i in 1:nrow(term)) {
            if (i == 1){
              clus_names <- term[i,4]
            }else{
              clus_names <- paste(clus_names, term[i,4],sep = ", ")
            }
          }
          #Create row that has term information, Trajectory Designation, how many 
          #clusters it is found in, and the names of those clusters
          term_row <- cbind(term[1,1:3],"Linear_Up",nrow(term),clus_names)
          colnames(term_row)[4:6] <- c("Direction","Num_Clusters","Clusters")
          #Add these to master list
          Divided_terms <- rbind(Divided_terms, term_row)
        }
        for (t in unique(lindown_unique$Gene.Set.ID)) {
          term <- lindown_unique %>% filter(Gene.Set.ID == t)
          for (i in 1:nrow(term)) {
            if (i == 1){
              clus_names <- term[i,4]
            }else{
              clus_names <- paste(clus_names, term[i,4],sep = ", ")
            }
          }
          term_row <- cbind(term[1,1:3],"Linear_Down",nrow(term),clus_names)
          colnames(term_row)[4:6] <- c("Direction","Num_Clusters","Clusters")
          Divided_terms <- rbind(Divided_terms, term_row)
        }
        for (t in unique(complex_unique$Gene.Set.ID)) {
          term <- complex_unique %>% filter(Gene.Set.ID == t)
          for (i in 1:nrow(term)) {
            if (i == 1){
              clus_names <- term[i,4]
            }else{
              clus_names <- paste(clus_names, term[i,4],sep = ", ")
            }
          }
          term_row <- cbind(term[1,1:3],"Complex",nrow(term),clus_names)
          colnames(term_row)[4:6] <- c("Direction","Num_Clusters","Clusters")
          Divided_terms <- rbind(Divided_terms, term_row)
        }
      }
    }
  }
  #Save master list once each enrichment term class has been cycled through
  write.csv(Divided_terms, "Aging Analysis Results/Enrichment/Divided_enrichment_w_Complex.csv")
}

#Grouping
gg <- "Gene_Group"
ont <- "Ontology"
path <- "Pathway"

#Subgroups
fb <- "Flybase"
slim2bp <- "SLIM2BP"
slim2mf <- "SLIM2MF"
slim2cc <- "SLIM2CC"
drsc <- "DRSC_GLAD"
rct <- "Reactome"

sorted_by_background(GROUPING = c(gg, ont, path),
                     SUBGROUPING = c(fb, slim2bp, slim2cc, slim2mf, drsc,rct))

sorted_unique <- read.csv("Aging Analysis Results/Enrichment/Divided_enrichment_w_Complex.csv", header = TRUE)
sorted_unique %>% dplyr::count(Direction)
#     Direction   n
# 1     Complex 125
# 2 Linear_Down 177
# 3   Linear_Up 293
sorted_unique %>% dplyr::count(Num_Clusters)
#   Num_Clusters   n
# 1            1 528
# 2            2  52
# 3            3  12
# 4            4   2
# 5            5   1

#So the majority of our unique terms are ones that are unique to only one cluster

#Divided Enrichment####
#Look at unique terms across the first divided in the dendrogram. This divided
#primarily is dividing based on if the trajectories primarily increasing or
#decreasing in expression.For some of the complex cluster's (Complex-1:5, 8, 10,11)
#this distinction is less clear, but for example comparing Complex 2 with Complex 11,
#rep curves, see complex 2 curve increases for most of lifespan then decreases, whereas
#complex 11 curve decreases for majority of the lifespan 

library(tidyverse)

#The code requires a specific directory setup.
# Master Directory/
#   deseq_results/
#     All_Analyses/
#       cluster_color_id_ref.csv
#   Aging Analysis Results/
#     Enrichment/
#       gene groups/
#         Complex.#/
#           Benjamini_Hochberg_enrich_genegroups_DRSC_GLAD_pangea_Complex.#.csv
#           Benjamini_Hochberg_enrich_genegroups_Flybase_pangea_Complex.#.csv
#       gene ontology/
#         Complex.#/
#           Benjamini_Hochberg_enrich_ontology_SLIM2BP_pangea_Complex.#.csv
#           Benjamini_Hochberg_enrich_ontology_SLIM2CC_pangea_Complex.#.csv
#           Benjamini_Hochberg_enrich_ontology_SLIM2MF_pangea_Complex.#.csv
#       pathways/
#         Complex.#/
#           Benhamini_Hochberg_enrich_pathway_Reactome_pangea_Complex.#.csv

#This code is very similar as Trajectory classification divide
#Only difference in towards end when seperating into Up and down
#groups vs LinearUp, LinearDown, Complex group

enrich_divide_background <- function(GROUPING, SUBGROUPING){
  #Read in cluster id information
  clust_desg <- read.csv("deseq_results/All Analyses/cluster_color_id_ref.csv", header = TRUE)
  #Separate clusters using Dendo_clus. The first divide break clusters into groups of 16 and 12.
  UP <- clust_desg %>% filter(Dendo_clus %in% 1:16)
  DOWN <- clust_desg %>% filter(!Dendo_clus %in% 1:16)
  #Create a master list to store divided terms
  Divided_terms <- NULL
  #For each enrichment term class
  for (g in GROUPING) {
    print(g)
    for (s in SUBGROUPING){
      print(s)
      #Want a list that has all the terms for that grouping and subgrouping
      enrich_terms <- NULL
      #For each cluster
      for (c in str_replace(clust_desg$Official_ID, "-", ".")) {
        if(g == "Gene_Group"){
          cluster_term_add <- paste("Aging Analysis Results/Enrichment/gene groups/",c,"/Benjamini_Hochberg_enrich_genegroups_",
                                    s,"_pangea_",c,".csv",sep = "")
        }else if (g == "Ontology"){
          cluster_term_add <- paste("Aging Analysis Results/Enrichment/gene ontology/",c,"/Benjamini_Hochberg_enrich_ontology_",
                                    s,"_pangea_",c,".csv",sep = "")
        }else if (g == "Pathway"){
          cluster_term_add <- paste("Aging Analysis Results/Enrichment/pathways/",c,"/Benjamini_Hochberg_enrich_pathway_",
                                    s,"_pangea_",c,".csv",sep = "")
        }
        if (file.exists(cluster_term_add)){
          #Read om Cluster's Significant PANGEA Terms
          cluster_term <- read.csv(cluster_term_add, header = TRUE)
          #Create a row with the term information and name of the cluster
          row <- cbind(cluster_term[,c(2:4)],str_replace(c,"\\.","-"))
          colnames(row)[4] <- c("Cluster")
          #bind the rows to the class specific list
          enrich_terms <- rbind(enrich_terms, row)
        }
      }
      #If there are terms in enrich_terms
      if(!is.null(enrich_terms)){
        #Divide the terms based on which side of the split they belong to
        up <- enrich_terms %>% filter(Cluster %in% UP$Official_ID)
        down <- enrich_terms %>% filter(Cluster %in% DOWN$Official_ID)
        #For each side of the split pull out only the terms that are unique to its side
        up_unique <- up %>% filter(!Gene.Set.ID %in% down$Gene.Set.ID)
        down_unique <- down %>% filter(!Gene.Set.ID %in% up$Gene.Set.ID)
        #Create a list of the clusters that have the unique terms.
        #Will do this for each trajectory so will only comment on Up
        for (t in unique(up_unique$Gene.Set.ID)) {
          term <- up_unique %>% filter(Gene.Set.ID == t)
          for (i in 1:nrow(term)) {
            if (i == 1){
              clus_names <- term[i,4]
            }else{
              clus_names <- paste(clus_names, term[i,4],sep = ", ")
            }
          }
          #Create a row that has term information, Side of Split, how many
          #clusters it is found, and names of those clusters
          term_row <- cbind(term[1,1:3],"Up",nrow(term),clus_names)
          colnames(term_row)[4:6] <- c("Direction","Num_Clusters","Clusters")
          #Add to master list
          Divided_terms <- rbind(Divided_terms, term_row)
        }
        for (t in unique(down_unique$Gene.Set.ID)) {
          term <- down_unique %>% filter(Gene.Set.ID == t)
          for (i in 1:nrow(term)) {
            if (i == 1){
              clus_names <- term[i,4]
            }else{
              clus_names <- paste(clus_names, term[i,4],sep = ", ")
            }
          }
          term_row <- cbind(term[1,1:3],"Down",nrow(term),clus_names)
          colnames(term_row)[4:6] <- c("Direction","Num_Clusters","Clusters")
          Divided_terms <- rbind(Divided_terms, term_row)
        }
      }
    }
  }
  #Save master list once each enrichment term class has been cycled through
  write.csv(Divided_terms, "Aging Analysis Results/Enrichment/Divided_enrichment.csv")
}

#Grouping
gg <- "Gene_Group"
ont <- "Ontology"
path <- "Pathway"

#Subgroups
fb <- "Flybase"
slim2bp <- "SLIM2BP"
slim2mf <- "SLIM2MF"
slim2cc <- "SLIM2CC"
drsc <- "DRSC_GLAD"
rct <- "Reactome"

enrich_divide_background(GROUPING = c(gg, ont, path), 
                         SUBGROUPING = c(fb, slim2bp, slim2cc, slim2mf, drsc,rct))

divide_unique <- read.csv("Aging Analysis Results/Enrichment/Divided_enrichment.csv", header = TRUE)
divide_unique %>% dplyr::count(Direction)
#   Direction   n
# 1      Down 269
# 2        Up 412
divide_unique %>% dplyr::count(Num_Clusters)
#   Num_Clusters   n
# 1            1 528
# 2            2 124
# 3            3  24
# 4            4   4
# 5            5   1

#Enrichment Term Count####
library(tidyverse)

#Create data table that has ever enriched gene with information of how 
#many clusters are enriched for it, and each enriched clusters
#has its own row. 

#The code requires a specific directory setup.
# Master Directory/
#   deseq_results/
#    All_Analyses/
#       cluster_color_id_ref.csv
#   Aging Analysis Results/
#     Enrichment/
#       gene groups/
#         enrich_genegroups_DRSC_GLAD.csv
#         enrich_genegroups_Flybase.csv
#       gene ontology/
#         enrich_SLIM2_BP.csv
#         enrich_SLIM2_CC.csv
#         enrich_SLIM2_MF.csv
#       pathways/
#         enrich_pathway_REACTOME.csv

term_count <- function(GROUPING, SUBGROUPING){
  term_count_df <- NULL
  for (g in GROUPING) {
    print(g)
    for (s in SUBGROUPING) {
      print(s)
      if(g == "Gene_Group"){
        pangea_add <- paste("Aging Analysis Results/Enrichment/gene groups/enrich_genegroups_",s,".csv",sep = "")
      }else if (g == "Ontology"){
        pangea_add <- paste("Aging Analysis Results/Enrichment/gene ontology/enrich_SLIM2_",s,".csv", sep="")
      }else if (g == "Pathway"){
        pangea_add <- paste("Aging Analysis Results/Enrichment/pathways/enrich_pathway_",s,".csv",sep = "")
      }
      if(file.exists(pangea_add)){
        pangea_term <- read.csv(pangea_add, header = TRUE)
        FDR_col <- pangea_term[substr(colnames(pangea_term),1,31)=="p.value..Benjamini...Hochberg.."]
        FDR_col <- FDR_col %>% replace(is.na(.),1)
        for (i in 1:nrow(pangea_term)) {
          term_sig <- NULL
          for (c in 1:ncol(FDR_col)) {
            print(c)
            if(FDR_col[i,c] < 0.05){
              term_sig <- append(term_sig, substring(colnames(FDR_col)[c],32))
            }
          }
          if(length(term_sig) != 0){
            submit_row <- cbind(remove_rownames(pangea_term[i,2:4]), length(term_sig), term_sig)
            colnames(submit_row) <- c("Gene.Set", "Gene.Set.ID","Gene.Set.Name", "Num_Clusters", "Clusters")
            term_count_df <- rbind(term_count_df, submit_row)
          }
        }
      }
    }
  }
  write.csv(term_count_df, "Aging Analysis Results/Enrichment/Term_Cluster_Count.csv")
}


#Grouping
gg <- "Gene_Group"
ont <- "Ontology"
path <- "Pathway"

#Subgroups
fb <- "Flybase"
slim2bp <- "BP"
slim2mf <- "MF"
slim2cc <- "CC"
drsc <- "DRSC_GLAD"
rct <- "Reactome"

term_count(GROUPING = c(gg, ont, path),
           SUBGROUPING = c(fb, slim2bp, slim2cc, slim2mf, drsc,rct))

#Ribosomal Protein Trajectories####
library(tidyverse)
library(rms)
library(RColorBrewer)
library(ggpubr)

#Identify genes that are associated with ribosomal related terms and examine
#their expression patterns.
#This code creates supplemental table M, and supplemental figures 13 and 14

#The code requires a specific directory setup.
# Master Directory/
#   deseq_results/
#     All Analyses/
#       Combo_Analysis_Clustered_Z_scores_k28.csv
#   Aging Analysis Results/
#     Enrichment/
#       gene groups/
#         enrich_genegroups_DRSC_GLAD.csv
#         enrich_genegroups_Flybase.csv
#       gene ontology/
#         enrich_SLIM2_BP.csv
#         enrich_SLIM2_CC.csv
#         enrich_SLIM2_MF.csv
#       pathways/
#         enrich_pathway_REACTOME.csv



#This function pulls out all genes associated either ribosomal terms that are listed
#in the TERMS variable. Genes may be pulled out multiple times. 
ribosomal_counts <- function(GROUPING, SUBGROUPING, TERMS){
  ribosomal_genes <- NULL
  for (g in GROUPING) {
    print(g)
    for (s in SUBGROUPING) {
      print(s)
      if(g == "Gene_Group"){
        pangea_add <- paste("Aging Analysis Results/Enrichment/gene groups/enrich_genegroups_",s,".csv", sep = "")
      }else if (g == "Ontology"){
        pangea_add <- paste("Aging Analysis Results/Enrichment/gene ontology/enrich_SLIM2_",s,".csv",sep = "")
      }else if (g == "Pathway"){
        pangea_add <- paste("Aging Analysis Results/Enrichment/pathways/enrich_pathway_",s,".csv", sep = "")
      }
      if (file.exists(pangea_add)){
        pangea <- read.csv(pangea_add, header = TRUE)
        pangea <- pangea[-c(1)]
        ribosomal_terms <- pangea %>% filter(Gene.Set.ID %in% TERMS)
        ribosomal_terms <- ribosomal_terms[c(1:3,228:255)]
        colnames(ribosomal_terms) <- gsub("species.specific.ID.","",colnames(ribosomal_terms))
        if (nrow(ribosomal_terms) != 0){
          for (r in 1:nrow(ribosomal_terms)) {
            #print(r)
            for (c in 4:ncol(ribosomal_terms)) {
              #print(c)
              genes <- str_split_1(ribosomal_terms[r,c],pattern = "; ")
              gene_set <- as.data.frame(cbind(ribosomal_terms[r,1:3], genes, colnames(ribosomal_terms)[c]))
              colnames(gene_set)[4:5] <- c("gene_id", "cluster")
              if(gene_set[1,4] != ""){
                ribosomal_genes <- rbind(ribosomal_genes, gene_set)
              }
            }
          }
        }
      }
    }
  }
  ribosomal_genes$Desg <- gsub("\\..*", "", ribosomal_genes$cluster)
  print(ribosomal_genes)
}

#Grouping
gg <- "Gene_Group"
ont <- "Ontology"
path <- "Pathway"

#Subgroups
fb <- "Flybase"
slim2bp <- "BP"
slim2mf <- "MF"
slim2cc <- "CC"
drsc <- "DRSC_GLAD"
rct <- "REACTOME"

GeneSetIDs <- c("GLAD:24600", 
                "GO:0042254", "GO:0003735","GO:0016072", 
                "R-DME-6791226", "R-DME-72312","FBgg0001149", "R-DME-73864","R-DME-72706","R-DME-72702","R-DME-8868773",
                "FBgg0000130","FBgg0000141", "FBgg0000200","FBgg0000071", "FBgg0000059", "FBgg0001149")

all_ribosomal_proteins <- ribosomal_counts(GROUPING = c(gg, ont, path), SUBGROUPING = c(fb, slim2bp, slim2cc, slim2mf, drsc, rct), 
                                           TERMS = GeneSetIDs)
head(all_ribosomal_proteins)
#             Gene.Set Gene.Set.ID                    Gene.Set.Name     gene_id      cluster       Desg
# 1 FlyBase Gene Group FBgg0000059 MITOCHONDRIAL RIBOSOMAL PROTEINS FBgn0035122 LinearDown.6 LinearDown
# 2 FlyBase Gene Group FBgg0000059 MITOCHONDRIAL RIBOSOMAL PROTEINS FBgn0035335 LinearDown.6 LinearDown
# 3 FlyBase Gene Group FBgg0000059 MITOCHONDRIAL RIBOSOMAL PROTEINS FBgn0038678 LinearDown.6 LinearDown
# 4 FlyBase Gene Group FBgg0000059 MITOCHONDRIAL RIBOSOMAL PROTEINS FBgn0011787 LinearDown.6 LinearDown
# 5 FlyBase Gene Group FBgg0000059 MITOCHONDRIAL RIBOSOMAL PROTEINS FBgn0037330 LinearDown.6 LinearDown
# 6 FlyBase Gene Group FBgg0000059 MITOCHONDRIAL RIBOSOMAL PROTEINS FBgn0036335 LinearDown.6 LinearDown

##Suplemental Table M####
#count the number of genes that are found for each of our terms. 
ribosomal_terms <- all_ribosomal_proteins %>% count(Gene.Set, Gene.Set.ID, Gene.Set.Name) %>% arrange(Gene.Set, Gene.Set.ID)
colnames(ribosomal_terms)[4] <- "Number of Genes"
head(ribosomal_terms)
#               Gene.Set Gene.Set.ID                        Gene.Set.Name Number of Genes
# 1 DRSC GLAD Gene Group  GLAD:24600                             Ribosome             123
# 2   FlyBase Gene Group FBgg0000059     MITOCHONDRIAL RIBOSOMAL PROTEINS              36
# 3   FlyBase Gene Group FBgg0000071 CYTOPLASMIC SMALL RIBOSOMAL PROTEINS              33
# 4   FlyBase Gene Group FBgg0000130                   RIBOSOMAL PROTEINS             118
# 5   FlyBase Gene Group FBgg0000141       CYTOPLASMIC RIBOSOMAL PROTEINS              80
# 6   FlyBase Gene Group FBgg0000200 CYTOPLASMIC LARGE RIBOSOMAL PROTEINS              46

write.csv(ribosomal_terms, "SupTableM_ribosomal_Terms.csv")

##Suplementary Figures 13 and 14####
#Create a table listing all genes associated with mitochondrial ribosomal proteins, and 
#a table that lists all other ribosomal related genes. 
mitochondrial_ribosomes <- all_ribosomal_proteins %>% filter(Gene.Set.ID == "FBgg0000059")
non_mitochondrial <- all_ribosomal_proteins %>% filter(Gene.Set.ID != "FBgg0000059") %>% filter(!gene_id %in% mitochondrial_ribosomes$gene_id)
non_mitochondrial <- unique(non_mitochondrial[4:6])

#Each table add a column designating if the gene is mitochondrial associated or not 
#and then combine both tables. 
non_mitochondrial$type <- "Not_Mitochondrial"
mitochondrial_sbuset <- cbind(mitochondrial_ribosomes[4:6],"Mitochondrial")
colnames(mitochondrial_sbuset)[4] <- "type"
ribosome_type <- rbind(non_mitochondrial, mitochondrial_sbuset)

#Read in Z scores of identified genes and combine with ribosomal type data frame
z_scores <- read.csv("deseq_results/All Analyses/Combo_Analysis_Clustered_Z_scores_k28.csv", header = TRUE)
z_scores_ribosomes <- left_join(ribosome_type, z_scores[-c(1,15,16)], by="gene_id")

z_scores_ribosomes$cluster <- gsub("\\.","-",z_scores_ribosomes$cluster)
z_scores_ribosomes$cluster <- factor(z_scores_ribosomes$cluster, levels = c("Complex-1","Complex-2","Complex-3","Complex-4",
                                                                            "Complex-5","Complex-6","Complex-7","Complex-8",
                                                                            "Complex-9", "Complex-10", "Complex-11","Complex-12",
                                                                            "Complex-13","Complex-14", "LinearUp-1","LinearUp-2",
                                                                            "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                                                            "LinearUp-7", "LinearUp-8", "LinearDown-1", "LinearDown-2",
                                                                            "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6"))

#Create list of cytoplasmic ribosomal proteins
cytoplasmic <- all_ribosomal_proteins %>% filter(Gene.Set.Name == "CYTOPLASMIC RIBOSOMAL PROTEINS")

#Mitochondrial vs Cytoplasmic 
#Plot expression trajecotires of only mitochondrial and cytoplasmic ribosomal proteins, regardless
#of cluster
mito_cyto_combine<- 
  z_scores_ribosomes %>% ggplot(aes(x=day, y=Z_score))+
  geom_smooth(data = z_scores_ribosomes %>% filter(gene_id %in% cytoplasmic$gene_id),method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE,aes(group=gene_id, color=type), linewidth=0.5)+
  geom_smooth(data = z_scores_ribosomes %>% filter(type!="Not_Mitochondrial"),method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE,aes(group=gene_id, color=type), linewidth=0.5)+
  scale_color_manual(values = brewer.pal(6,"Paired")[c(4,6)],
                     labels = c("FBgg0000059: Mitochondrial Ribosomal Proteins",
                                "FBgg0000141: Cytoplasmic Ribosomal Proteins"),
                     name="Ribosomal Protein")+
  ylab("Z score")+
  xlab(("Day"))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(text = element_text(size=12))

#Plot expression trajectories of only mitochondrial and cytoplasmic ribosomal proteins,
#and facet based on cluster
mito_cyto_clus <- 
  z_scores_ribosomes %>%filter(gene_id %in% cytoplasmic$gene_id |type!="Not_Mitochondrial")%>% ggplot(aes(x=day, y=Z_score))+
  facet_wrap(~cluster)+
  geom_smooth(data = z_scores_ribosomes %>% filter(gene_id %in% cytoplasmic$gene_id),method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE,aes(group=gene_id, color=type), linewidth=0.5)+
  geom_smooth(data = z_scores_ribosomes %>% filter(type!="Not_Mitochondrial"),method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE,aes(group=gene_id, color=type), linewidth=0.5)+
  scale_color_manual(values = brewer.pal(6,"Paired")[c(4,6)],
                     labels = c("FBgg0000059: Mitochondrial Ribosomal Proteins",
                                "FBgg0000141: Cytoplasmic Ribosomal Proteins"),
                     name="Ribosomal Protein")+
  ylab("Z score")+
  xlab(("Day"))+
  theme_bw()+
  theme(text = element_text(size=12))

#merge both figures and save
ribosomes_cyto_mito<-ggarrange(mito_cyto_combine, mito_cyto_clus, ncol=2, common.legend = TRUE, legend = "bottom")
ggsave("Plots/SupFig13_Mito_Cyto.png",ribosomes_cyto_mito,width = 8.5, height = 5, units = "in")

#Plot expression trajectories of mitochondrial proteins and all other ribosomal associated proteins. 
mito_all_combine<- 
  z_scores_ribosomes %>% ggplot(aes(x=day, y=Z_score))+
  geom_smooth(data = z_scores_ribosomes ,method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE,aes(group=gene_id, color=type), linewidth=0.5)+
  geom_smooth(data = z_scores_ribosomes %>% filter(type!="Not_Mitochondrial"),method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE,aes(group=gene_id, color=type), linewidth=0.5)+
  scale_color_manual(values = brewer.pal(6,"Paired")[c(4,5)],
                     labels = c("FBgg0000059: Mitochondrial Ribosomal Proteins",
                                "All Other Ribosomal Related Genes"),
                     name="Ribosomal Protein")+
  ylab("Z score")+
  xlab(("Day"))+
  theme_bw()+
  theme(text = element_text(size=12))

#Plot expression trajectories of mitochondrial proteins and all other ribosomal associated proteins,
#and facet by cluster. 
mito_all_clus <- 
  z_scores_ribosomes %>% ggplot(aes(x=day, y=Z_score))+
  facet_wrap(~cluster, ncol=6)+
  geom_smooth(method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE, show.legend = TRUE,aes(group=gene_id,color=type), linewidth=0.5)+
  scale_color_manual(values = brewer.pal(6,"Paired")[c(4,5)],
                     labels = c("FBgg0000059: Mitochondrial Ribosomal Proteins",
                                "All Other Ribosomal Related Genes"),
                     name="Ribosomal Protein")+
  ylab("Z score")+
  xlab(("Day"))+
  theme_bw()+
  theme(text = element_text(size=8.5),
        legend.position = "none")

#merge both figures and save
ribosomes_all_mito<-ggarrange(mito_all_combine, mito_all_clus, ncol=2, common.legend = TRUE, legend = "bottom")
ggsave("Plots/SupFig14_Mito_all.png",ribosomes_all_mito,width = 9, height = 5, units = "in")
