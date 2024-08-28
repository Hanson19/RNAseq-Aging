#Introduction####
#Title: Tximeta transcript quantification import
#Packages:
#install.packages("BiocManager")
#BiocManager::install("tximeta")
library(tximeta)
suppressPackageStartupMessages(library(SummarizedExperiment))

#Read in sample table
sample_table <- read.csv("sample_table_directory_RNAseq_Aging.csv", header = T,
                         colClasses = c("character","character","numeric","numeric",
                                        "factor"))
#Read in a sample table that has for each sample, the file address to its quant.sf files,
#columns with the samples name, the day it was collected, the population survivorship, and what replicate.
colnames(sample_table)
# [1] "files"        "names"        "day"          "survivorship" "replicate"

se <- tximeta(sample_table)

# Look at the summarized experiment generated####
colData(se) # this is the sample table
assayNames(se) # has counts/abundance/length and bootstrap-based inferential replicates
rowRanges(se) # transcript detail
seqinfo(se) # detail on chromosomes


# Add info about exons per transcript
se.exons <- addExons(se)
rowRanges(se.exons)[[1]] # Look at exon info


# Summarize to gene level
# (Get following warning when run this:
#  call dbDisconnect() when finished working with a connection
#  According to https://support.bioconductor.org/p/9137815/
#  this is ignorable)
gse <- summarizeToGene(se)
rowRanges(gse) # look at object
names(metadata(gse)) # examine metadata
str(metadata(gse)[["quantInfo"]]) # examine structure of any metadata object


# Save all objects as RDS files
#  Note these files are big (se and se_ex are ~1Gb,
#  and gse is 350Mb) and slow to generate
saveRDS(se,file="se_salmon_tximeta.RDS")
saveRDS(se.exons,file="se_ex_salmon_tximeta.RDS")
saveRDS(gse,file="gse_salmon_tximeta.RDS")