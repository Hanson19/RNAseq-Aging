# RNAseq-Aging

#### Scripts Available:

* Identifying_and_Clustering.R
  * Identifying genes whose expression changes with aging and clustering those
  genes using their expression trajectories.
  * Determining if a cluster has a complex or linear trajectory
  * Creating Figures 1-2, Supplementary Figures 3 and 4, and Supplementary tables D and E
* Designation Analysis ID Enrichment.R
  * Comparing Clusters and trajectories gene analysis identification.
  * Creating Supplementary Figures 5 and 6
  * Does require files made in Identifying_and_Clustering.R
* Cluster Enrichment Analysis.R
  * Analysis of PANGEA results and creation of summary tables
  * Does require files made in Identifying_and_Clystering.R
* Young v Old Analysis.R
  * Young v Old Analysis using samples from Days 3 and 6 as the young, and Day 59's samples the old.
  * Compared identified genes in analysis with those identified in the multi-timepoint analysis. 
  * Looked at the distribution of YvO direction calls across clusters and trajectory designations.
  * Creating Figure 4 and Supplementary Figure 8
  * Requires files made in Identifying_and_Clustering.R
* Young v Alternative Old Analysis.R
  * Young v Old Analyses using alternative old samples which come from two consecutive collection days.
  * Compared gene identification and direction calls of the alternative analyses and the original
  Day 59 analysis as well as the multi-timepoint analysis.
  * Creating Supplementary Figures 9 and 10
  * Requires files made in Identifying_and_Clustering.R
* Comparison with Published Aging Studies.R
  * Comparing our identified multi-timepoint genes with 9 data frames from previously published Young v Old aging expression analyses. 
  * Creating Supplementary Figures 11 and 12, and Supplementary Table L
* Misc Supplemental Figure Code.R
  * Creating Supplementary Figures 2 and 7
  * Require files made in Identifying_and_Clustering.R
* Aging_Survivorship_Code.R
  * Kaplan-Meier Survivorship curve
  * Creating Supplementary Figure 1
  
