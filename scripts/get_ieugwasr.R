#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
get_ieugwasr.R
===============

Author: Antonio Berlanga
Release: |version|
Date: |today|


Purpose
=======

Extract GWAS summary statistics for outcome data for two sample Mendelian randomisation using the R package
ieugwasr.


Usage and options
=================

Usage: get_ieugwasr.R (--instruments <file>) (--pmid <pmid>) [options]
       get_ieugwasr.R [-h | --help]

Options:
  --instruments <file>   Single column file without header with instruments to search for in outcome data in IEU GWAS database
  --pmid <pmid>          PubMed ID of study to search for as outcome data
  -O <OUTPUT_FILE>       Output file prefix
  --session              R session if to be saved [default: FALSE]
  -h --help              Show this screen

Input:

   Single column file without a header with SNP instruments to search for in outcome data. PMID of the outcome file is needed to specify which study to search against.

   Additionally, an access token to the IEU GWAS database need to be provided. This is a manual step that needs to be done in an interactive session. If batch running generate the token and copy across the "ieugwasr_oauth" folder to the working directory. For more information see:
   https://mrcieu.github.io/ieugwasr/articles/guide.html#get-list-of-a-specific-study

Output:

    If the PMID is in theIEU GWAS database a tab separated file with headers with SNPs matching instruments.

Command example:

  get_ieugwasr.R --instruments my_file.txt --pmid 26502338 -O SLE_26502338

Requirements:

	episcout
	ieugwasr

Documentation
=============

    For more information see:

    |url|

' -> doc

# Load docopt:
library(docopt, quietly = TRUE)
# Retrieve the command-line arguments:
args <- docopt(doc)

# Print args given to screen:
str(args)

# options(echo = TRUE) # print commands for debugging
######################

######################
# This function allows other R scripts to obtain the path to a script directory
# (ie where this script lives). Useful when using source('some_script.R')
# without having to pre-specify the location of where that script is.
# This is taken directly from:
# How to source another_file.R from within your R script molgenis/molgenis-pipelines Wiki
# https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
# Couldn't find a licence at the time (12 June 2018)
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }

  if (!is.null(this.file)) return(dirname(this.file))

  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from = 1, length.out = length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))

  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}
Rscripts_dir <- LocationOfThisScript()
print('Location where this script lives:')
Rscripts_dir
# R scripts sourced with source() have to be in the same directory as this one
# (or the path constructed appropriately with file.path) eg:
#source(file.path(Rscripts_dir, 'moveme.R')) #, chdir = TRUE)
######################

######################
# Import libraries
library(episcout)
library(ieugwasr)
# source functions from a different R script:
# source(file.path(Rscripts_dir, 'ggtheme.R'))
######################

######################
# Required input files

##########
# Read instruments data:
if (!is.null(args[['--instruments']])) { # for docopt this will be NULL or chr, if boolean
  # remove is.null function and test with ==
  instruments_file <- as.character(args[['--instruments']])
  instruments <- epi_read(instruments_file, header = FALSE)
  instruments <- instruments$V1
  # For tests:
  # setwd('~/Documents/quickstart_projects/projects/MR_pipeline_runs/pipeline_MR_tests/ieugwasr_tests')
  # instruments_file <- 'cis_pqtl_instruments_SNPs_only.txt'
  # instruments <- epi_read(instruments_file, header = FALSE)
  # instruments <- instruments$V1
  } else {
  # Stop if arguments not given:
  print('You need to provide an instruments input file.')
  stopifnot(!is.null(args[['--instruments']]))
}

print('File being used for instruments data: ')
print(instruments_file)
##########

##########
# Read PMID argument:
if (!is.null(args[['--pmid']])) { # for docopt this will be NULL or chr, if boolean
  # remove is.null function and test with ==
  pmid <- as.character(args[['--pmid']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/projects/MR_pipeline_runs/pipeline_MR_tests/ieugwasr_tests')
  # pmid <- '26502338'
  } else {
  # Stop if arguments not given:
  print('You need to provide a Pubmed ID to search for.')
  stopifnot(!is.null(args[['--pmid']]))
  }

print('PMID passed: ')
print(pmid)
##########

##########
# Set output file names:
if (is.null(args[['-O']])) { # arg is NULL
  # Use PMID as default:
  output_prefix <- sprintf('%s_in_%s', instruments_file, pmid)
  print('Output file prefix not given. Using: ')
  print(output_prefix)
  } else {
    output_prefix <- as.character(args[['-O']])
  # output_file_prefix <- 'testing'
  print(sprintf('Output file name provided: %s', output_prefix))
  print(output_prefix)
  }
##########
######################

######################
# Sign in, check status:
# ieugwasr::get_access_token()

if (is.null(ieugwasr::check_access_token())) {
  stop('An access token to the IEU database was not found. This needs to be done manually from the working directory (or copy across). From an R session run:
  ieugwasr::get_access_token()
  For more information see:
        https://mrcieu.github.io/ieugwasr/articles/guide.html#get-list-of-a-specific-study')
  }
# api_status()
######################

######################
##########
# Check what studies are available:
data_available <- ieugwasr::gwasinfo()

# Search outcome by PMID:
# pmid <- '31604244'
get_study <- which(data_available$pmid == pmid)
if (length(get_study) == 0) {
  stop('PMID not found in IEU GWAS database. Exiting.')
  }
# as.data.frame(data_available[get_study, ])
outcome_info <- data_available[get_study, ]
# as.data.frame(outcome_info)
ieu_id <- outcome_info$id
# ieu_id

# Get outcome variants:
# This will get matching SNPs, check for proxies, align alleles, etc.:
# TO DO: pass as args
get_outcome_SNPs <- ieugwasr::associations(variants = instruments,
                                           id = ieu_id,
                                           proxies = 1,
                                           r2 = 0.8,
                                           align_alleles = 1,
                                           palindromes = 1,
                                           maf_threshold = 0.01,
                                           access_token = check_access_token()
                                           )
# get_outcome_SNPs
##########

##########
# Rename columns so that they match TwoSampleMR package:
# TwoSampleMR columns:
# "SNP",
# "se",
# "pval",
# "beta",
# "samplesize",
# "effect_allele",
# "other_allele",
# "eaf",
# "Phenotype",

# Not in ieugwasr output:
# "units",
# "ncase",
# "ncontrol",
# "gene",
# "id"

# colnames(get_outcome_SNPs)
colnames(get_outcome_SNPs)[which(colnames(get_outcome_SNPs) == 'rsid')] <- 'SNP'
colnames(get_outcome_SNPs)[which(colnames(get_outcome_SNPs) == 'ea')] <- 'effect_allele'
colnames(get_outcome_SNPs)[which(colnames(get_outcome_SNPs) == 'nea')] <- 'other_allele'
colnames(get_outcome_SNPs)[which(colnames(get_outcome_SNPs) == 'p')] <- 'pval'
colnames(get_outcome_SNPs)[which(colnames(get_outcome_SNPs) == 'n')] <- 'samplesize'
colnames(get_outcome_SNPs)[which(colnames(get_outcome_SNPs) == 'trait')] <- 'Phenotype'
##########

##########
# Write to disk:
output_file <- sprintf('%s.tsv', output_prefix)
print('Saving results to:')
print(output_file)
epi_write(get_outcome_SNPs, output_file)
##########
######################

######################
# End IEUGWASR session:
# ieugwasr::revoke_access_token()
######################

######################
# The end
# Filename to save current R session, data and objects at the end:
if (isTRUE(args[['--session']])) {
  save_session <- sprintf('ieugwasr_%s.RData', output_prefix)
  print(sprintf('Saving an R session image as: %s', save_session))
  save.image(file = save_session, compress = 'gzip')
} else {
  print('Not saving an R session image, this is the default. Specify the --session option otherwise')
}

print('Finished successfully.')
sessionInfo()
q()
# Next: run the script for xxx
######################
