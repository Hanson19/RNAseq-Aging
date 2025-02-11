# RNAseq-Aging

Data and Code for multi-timepoint RNAseq analysis. 

The R and shell scripts used for the multi-timepoint RNAseq, enrichment, and young vs old analyses are provided in the scripts directory. 
All R scripts (exception Tximeta transcript quantification import.R) can be run with the provided data in data_frames directory.
These scripts begin with an introduction that includes a brief description of the scripts function and are divided into sections which begin 
with what packages are required and the expected directory structure/files. 

Required data frames to run these scripts are available in data_frames directory, and 
are separated into different directories. These directories have the same name as those
expected in the directory structure. 

Below are listed and described the available scripts and data_frame directories. 

#### Scripts Available:

* Identifying_and_Clustering.R
  * Identifying genes whose expression changes with aging and clustering those
  genes using their expression trajectories.
  * Determining if a cluster has a complex or linear trajectory
  * Creating Figures 2-4, Supplementary Figures 2 and 3, and Supplementary tables D and E
* Designation Analysis ID Enrichment.R
  * Comparing Clusters and trajectories gene analysis identification.
  * Creating Supplementary Figures 4 and 5
  * Does require files made in Identifying_and_Clustering.R
* Cluster Enrichment Analysis.R
  * Analysis of PANGEA results and creation of summary tables
  * Creating Figure 5, and Supplementary table L
  * Does require files made in Identifying_and_Clystering.R
* Young v Old Analysis.R
  * Young v Old Analysis using samples from Days 3 and 6 as the young, and Day 59's samples the old.
  * Compared identified genes in analysis with those identified in the multi-timepoint analysis. 
  * Looked at the distribution of YvO direction calls across clusters and trajectory designations.
  * Creating Figure 6 and Supplementary Figure 7
  * Requires files made in Identifying_and_Clustering.R
* Young v Alternative Old Analysis.R
  * Young v Old Analyses using alternative old samples which come from two consecutive collection days.
  * Compared gene identification and direction calls of the alternative analyses and the original
  Day 59 analysis as well as the multi-timepoint analysis.
  * Creating Supplementary Figures 8 and 9
  * Requires files made in Identifying_and_Clustering.R
* Comparison with Published Aging Studies.R
  * Comparing our identified multi-timepoint genes with 9 data frames from previously published Young v Old aging expression analyses. 
  * Creating Supplementary Figures 10 and 11, and Supplementary Table M
* Misc Supplemental Figure Code.R
  * Creating Supplementary Figures 1 and 6
  * Require files made in Identifying_and_Clustering.R
* Aging_Survivorship_Code.R
  * Kaplan-Meier Survivorship curve
  * Figure 1
* Tximeta transcript quantification import.R
  * Script we used to get from quant.sf files for individual samples and
  summarize the results to a gene level.
  * Cannot run without sample table and the sample quant files.
* aging_salmon_alignment.sh
  * Salmon pseudo alignment shell script 
  
#### Data frame directories:

* Comparison with Other Studies
  * Have available for each of the studies the validated Gene ID list from flybase and
  a file that allowed us to get the studies significant genes and direction calls.
  * Also available validated Gene IDs for significant multi-timepoint analysis genes. 
* Old v Young sim
  * Has the gse_salmon_tximeta file for the OvY analysis using Day 59 as old samples
  * Has an additional directory "Alternative old Code", which has these files for the 
  alternative old analyses
* Quant Files
  * gse_salmon_tximeta file for multi-timepoint anlaysis.
  * This is the output from the tximeta transcript quantification import.R script.
* Survivorship 
  * Death matrix needed for Kaplan-Meier Survivorship curve
  * Output of Aging_Survivorship_Code which is needed for Identifying_and_Clustering.R
* gene groups
  * PANGEA results tables for DRSC GLAD and Flybase gene groups
* gene ontology
  * PANGEA results tables for SLIM2 BP, MF, and CC ontology terms
* pathways
  * PANGEA results table for REACTOME pathways
