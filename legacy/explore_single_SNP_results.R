######################
# Explore results from single SNP analysis
# April 2020
# Antonio
######################


######################
setwd('/Users/antoniob/Documents/quickstart_projects/projects/cytokines_blood_traits_MR/results/parsa_30_Apr_2020/results_summary/')
######################

######################
library(episcout)
library(data.table)
######################

######################
# Read files:
single_SNP_results <- epi_read('single_SNP.summary_tsv')
epi_head_and_tail(single_SNP_results)
colnames(single_SNP_results)
dim(single_SNP_results)
str(single_SNP_results)

# Number of unique SNPs:
length(unique(single_SNP_results$SNP))

# Unique number of exposure - outcomes pairs:


# Number of significant results 0.05:
length(which(single_SNP_results$p < 0.05))
summary(single_SNP_results$p)
summary(single_SNP_results$b)

#

######################





######################
# Write files to disk:
# Abbreviations in this file will have duplicates, use filenames for now:
epi_write()
epi_write()
######################

######################
# The end
sessionInfo()
q()
######################
