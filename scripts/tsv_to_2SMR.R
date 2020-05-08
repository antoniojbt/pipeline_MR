##########
# pipeline_MR utility scripts
# April 2020
# Antonio
# Format GWAS summary data for TwoSampleMR package
##########

##########
# setwd('/Users/antoniob/Documents/github.dir/AntonioJBT/pipeMR/tests/testthat/')
##########

##########
# Set up command-line options
options(echo = TRUE) # see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# input_file <- as.character(args[1])
##########

##########
library(episcout)
##########

##########
# Read in full file
# input_file <- ''
input_file <- as.character(args[1])
# Header false if using grep'd files other T
mr_raw_data <- episcout::epi_read(input_file, header = F)
epi_head_and_tail(mr_raw_data)
epi_head_and_tail(mr_raw_data, last_cols = T)
# epi_clean_get_dups(mr_raw_data)
# unique(mr_raw_data)
str(mr_raw_data)
# View(mr_raw_data)
##########

##########
# TO DO: skip if full file
# Get column names and add back as lost during grepping:
colnames(mr_raw_data)
colnames_to_change <- c()

colnames(mr_raw_data) <- colnames_to_change
colnames(mr_raw_data)
epi_head_and_tail(mr_raw_data)
# View(mr_raw_data)
##########

##########
# Match names to TwoSampleMR, keep exp and out as separate files until
# harmonising

# TwoSampleMR columns:
# "Phenotype",
# "SNP",
# "beta",
# "se",
# "eaf",
# "effect_allele",
# "other_allele",
# "pval",

# Missing:
# "samplesize",
# "units",
# "ncase",
# "ncontrol",
# "gene",
# "id"

colnames(mr_raw_data)
colnames(mr_raw_data)[which(colnames(mr_raw_data) == 'ID')] <- 'SNP'
colnames(mr_raw_data)[which(colnames(mr_raw_data) == 'EFFECT')] <- 'beta'
colnames(mr_raw_data)[which(colnames(mr_raw_data) == 'SE')] <- 'se'
colnames(mr_raw_data)[which(colnames(mr_raw_data) == 'REF')] <- 'effect_allele'
colnames(mr_raw_data)[which(colnames(mr_raw_data) == 'ALT')] <- 'other_allele'
colnames(mr_raw_data)[which(colnames(mr_raw_data) == 'P')] <- 'pval'
colnames(mr_raw_data)[which(colnames(mr_raw_data) == 'N')] <- 'samplesize'
colnames(mr_raw_data)
##########


##########
# Get phenotype from file names and add as column for each SNP:
# ideally only split at '.'
pheno_name <- strsplit(input_file, split = '[.]')[[1]][1]
§pheno_name

mr_raw_data$Phenotype <- pheno_name
mr_raw_data[, 'Phenotype']

# Get reference allele frequency from ALT_FREQ:
mr_raw_data$eaf <- 1 - mr_raw_data$ALT_FREQ
mr_raw_data[, c('ALT_FREQ', 'eaf')]
##########


##########
# Write single file to disk in TwoSampleMR format
# If running directly:
# outfile_name <- epi_output_name(input_name = input_file, suffix = '.2SMR.tsv')
# If running from cgatcore pipeline:
outfile_name <- as.character(args[2])
epi_write(mr_raw_data, file_name = outfile_name)
##########

##########
sessionInfo()
q()
##########
