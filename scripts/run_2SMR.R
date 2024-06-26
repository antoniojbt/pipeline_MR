#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
run_2SMR
===============

Author: Antonio Berlanga
Release: |version|
Date: |today|


Purpose
=======

Perform analysis for two sample Mendelian randomisation using the R package
TwoSampleMR. Analysis for each SNP individually is ran by default, further tests can be added in options.


Usage and options
=================

Usage: run_2SMR.R (--exposure <exposure_file>) (--outcome <outcome_file>) [options]
       run_2SMR.R [-h | --help]

Options:
  --exposure <exposure_file>      Exposure input file name
  --outcome <outcome_file>        Outcome input file name
  -O <OUTPUT_FILE>                Output file prefix for exposure ("_on_" is added so that eg bmi_on_CHD is used as a prefix)
  --mr-methods <chr>              MR methods "none", "main" or "all" [default: main]
  --run-radial                    Run MR IVW radial method [default: FALSE]
  --run-raps                      Run MR RAPS method [default: FALSE]
  --stop-clumping                 Stop clumping of SNPs [default: FALSE]
  --local-LD                      Run LD operations locally. Requires bed, bim, fam files [default: false]
  --clump-kb <num>                TwoSampleMR clump_data() parameter [default: 10000]
  --clump-r2 <num>                TwoSampleMR clump_data() parameter [default: 0.001]
  --clump-p1 <num>                TwoSampleMR clump_data() parameter [default: 1]
  --clump-p2 <num>                TwoSampleMR clump_data() parameter [default: 1]
  --action <int>                  TwoSampleMR action to take for harmonising [default: 2]
  --egger-robust <BOOL>           TwoSampleMR mr_egger() parameter [default: FALSE]
  --egger-penalized <BOOL>        TwoSampleMR mr_egger() parameter [default: FALSE]
  --egger-distribution <chr>      TwoSampleMR mr_egger() parameter [default: normal]
  --egger-alpha <num>             TwoSampleMR mr_egger() parameter [default: 0.05]
  --OUTLIERtest <BOOL>            TwoSampleMR mr_egger() parameter [default: TRUE]
  --DISTORTIONtest <BOOL>         TwoSampleMR mr_egger() parameter [default: TRUE]
  --NbDistribution <int>          TwoSampleMR mr_egger() parameter [default: 1000]
  --SignifThreshold <num>         TwoSampleMR mr_egger() parameter [default: 0.05]
  --outcome-type <chr>            For Steiger test, "quant" or "binary" [default: binary]
  --session                       R session if to be saved [default: FALSE]
  -h --help                       Show this screen

Input:

    A tab separated file with headers which conform to TwoSampleMR
    requirements.

Output:

Various plots and tables for 2SMR.

Notes:

To run LD locally a "plink" executable must be in the PATH as well as ".bed", ".bim" and ".fam" files in the
working directory.

Command example:

  run_2SMR.R --exposure exposure_bmi_TwoSampleMR.tsv --outcome outcome_chd_TwoSampleMR.tsv --run-radial -O test --egger-robust TRUE --egger-penalized TRUE --stop-clumping --OUTLIERtest FALSE --session --mr-methods all

Output files when using "-O test" will start as test_on_outcome_chd_TwoSampleMR

Requirements:

	episcout
	dplyr
	ggplot2
	cowplot
	psych
	RadialMR
	mr.raps
	MRPRESSO
	MendelianRandomization
	TwoSampleMR
    ieugwasr

Documentation
=============

    Defaults are those set by the relevant packages called (TwoSampleMR, MendelianRandomization, etc.)
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
# TO DO:

# convert to functions and move them to pipeMR

# function to make table from:
# mr_single_SNP_pqtl_xxx.txt
# mr_results_pqtl_xxx.txt
# mr_heterogeneity_pqtl_xxx.txt
# mr_egger_i2_pqtl_xxx.txt
# mr_loo_IVW_pqtl_xxx.txt
# mr_presso_pqtl_xxx.txt
# mr_pleiotropy_pqtl_xxx.txt
# mr_steiger_pqtl_xxx.txt

# Other results are:
# mr_plots_pqtl_xxx.svg

# pass as arg:
# TwoSampleMR::mr_method_list() # hardcoded, currently all but radial and raps
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
message('Location where this script lives:')
Rscripts_dir
# R scripts sourced with source() have to be in the same directory as this one
# (or the path constructed appropriately with file.path) eg:
#source(file.path(Rscripts_dir, 'moveme.R')) #, chdir = TRUE)
######################

######################
# Import libraries
library(episcout)
library(dplyr)
library(ggplot2)
library(cowplot)

library(psych)
library(RadialMR)
library(mr.raps)
library(MRPRESSO)

library(MendelianRandomization)
library(TwoSampleMR)
library(ieugwasr)

# source functions from a different R script:
# source(file.path(Rscripts_dir, 'ggtheme.R'))
######################

######################
# Required input files

##########
# Read exposure data:
if (!is.null(args[['--exposure']])) { # for docopt this will be NULL or chr, if boolean
	                            # remove is.null function and test with ==
  exposure_file <- as.character(args[['--exposure']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/projects/MR_pipeline_runs/pipeline_MR_tests')
  # setwd('~/Documents/quickstart_projects/projects/cytokines_blood_traits_MR/results/ville_pQTL_on_AID/')
  # exposure_file <- 'exposure_bmi_TwoSampleMR.tsv'
  # exposure_file <- 'cis_pqtl_instruments.2SMR.tsv'
  exp_data <- TwoSampleMR::read_exposure_data(filename = exposure_file,
                                              sep = '\t',
                                              clump = FALSE,
                                              phenotype_col = 'Phenotype'
                                              )
  # colnames(exp_data)
  # head(exp_data)
  # dim(exp_data)

} else {
  # Stop if arguments not given:
  warning('You need to provide an exposure input file. This has to be tab separated with headers in the format for TwoSampleMR.')
  stopifnot(!is.null(args[['--exposure']]))
}

message('File being used for exposure data: ')
message(exposure_file)
##########

##########
# Read outcome data:
if (!is.null(args[['--outcome']])) { # for docopt this will be NULL or chr, if boolean
  # remove is.null function and test with ==
  outcome_file <- as.character(args[['--outcome']])
  # For tests:
  # outcome_file <- 'outcome_chd_TwoSampleMR.tsv'
  # outcome_file <- 'outcome_data_available/AS_23749187.tsv'

  out_data <- TwoSampleMR::read_outcome_data(filename = outcome_file,
                                             sep = '\t',
                                             phenotype_col = 'Phenotype'
                                             )
  # colnames(out_data)
  # head(out_data)
  # dim(out_data)
  } else {
  # Stop if arguments not given:
  warning('You need to provide an outcome input file. This has to be tab separated with headers in the format for TwoSampleMR.')
  stopifnot(!is.null(args[['--outcome']]))
}
message('File being used for outcome data: ')
message(outcome_file)
##########

##########
# Set output file names:
if (is.null(args[['-O']])) { # arg is NULL
  # Split infile names at the last '.':
  outfile_name_exp <- strsplit(exposure_file, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  outfile_name_out <- strsplit(outcome_file, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  output_file_prefix <- smessagef('%s_on_%s', outfile_name_exp, outfile_name_out)
  message('Output file prefix not given. Using: ')
  message(output_file_prefix)
} else {
  outfile_name_exp <- as.character(args[['-O']])
  # output_file_prefix <- 'testing'
  outfile_name_out <- strsplit(outcome_file, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  output_file_prefix <- sprintf('%s_on_%s', outfile_name_exp, outfile_name_out)
  message('Output file prefix provided, using: ')
  message(output_file_prefix)
}
##########
######################

######################
# Clump and harmonise

##########
# Clump SNPs:
# Set up args:
clump_kb <- as.numeric(args[['--clump-kb']])
clump_r2 <- as.numeric(args[['--clump-r2']])
clump_p1 <- as.numeric(args[['--clump-p1']])
clump_p2 <- as.numeric(args[['--clump-p2']])
# clump_kb <- 10000
# clump_r2 <- 0.001
# clump_p1 <- 1
# clump_p2 <- 1

if (isTRUE(args[['--stop-clumping']])) {
  message('--stop-clumping option used, skipping clumping.')
  } else if (isTRUE(args[['--local-LD']])) {
      message('Running LD locally')
      plink_files <- list(list.files(pattern = "\\.bed$"),
                          list.files(pattern = "\\.bim$"),
                          list.files(pattern = "\\.fam$")
                          )
      if (length(plink_files) != 3) {
          stop('Plink files not present, bim, bed and fam files required to run LD locally. See:
                   https://mrcieu.github.io/ieugwasr/articles/local_ld.html'
                   ) } else {
                         message('Clumping SNPs using local reference')
                         exp_data <- ieugwasr::ld_clump(dplyr::tibble(rsid = exp_data$SNP,
                                                                      pval = exp_data$pval,
                                                                      id = exp_data$Phenotype,
                                                                      clump_kb = clump_kb,
                                                                      clump_r2 = clump_r2,
                                                                      clump_p = clump_p1
                                                                      ),
                                                                      plink_bin = 'plink',
                                                                      bfile = '.'
                                                                      )
                   }
  } else {
    message('Clumping SNPs through MRCIEU API')
    exp_data <- TwoSampleMR::clump_data(dat = exp_data,
                                        clump_kb = clump_kb,
                                        clump_r2 = clump_r2,
                                        clump_p1 = clump_p1,
                                        clump_p2 = clump_p2
                                        )
    }
# head(exp_data)
# dim(exp_data)

# TO DO: get LD matrix through API:
# ieugwasr::ld_matrix(b$variant)
# or locally:
# ieugwasr::ld_matrix(exp_data$SNP,
#                     plink_bin = 'plink',
#                     bfile = '.'
#                     )

# TO DO: get LD proxies:
# See:
# https://mrcieu.github.io/gwasvcf/
# https://mrcieu.github.io/gwasvcf/articles/guide.html#ld-proxies-1
##########

##########
# TO DO:
# Calculate F-statistic for each instrument:
# For each of the instruments, if you have the SNP-exposure data then
# qf(pval, df1=1, df2=samplesize-1, lower.tail=FALSE)
# Also can be:
# We tested the statistical significance of the association between the instrument and CRP with the following formula:
# F=R2(n−1−k)(1−R2)×k 
#  We estimated the variance explained in serum amounts of CRP by using the formula:
# (2 × MAF(1 − MAF)β2)/var(CRP), where β is the estimated effect of the individual variants on
# CRP23 and var(CRP) is the estimated variance in natural-log-transformed CRP in the RS-I cohort.

# Save table with basic info for each instrument
# ######
# - Script: Test instrument strength with partial F statistics and r2 variance for each SNP, keep those > 10. This will exclude weak instruments and reduce weak instrument bias. A separate power calculation is needed however to test specific hypotheses.
#     + Approximate F-statistic calculation
#     + Exact F-statistic calculation
#     + Calculate pooled F-statistic
# eg:
# For each of the instruments, if you have the SNP-exposure data then
# qf(pval, df1=1, df2=samplesize-1, lower.tail=FALSE)


# # From Xiyun:
# # Variance explained and F-statistics
# R2<-rep(0,275)
# f_stats<-rep(0,275)
# for(i in 1:275){
#   R2[i] <-(2*(clumppped_data[i, 5])^2*clumppped_data[i, 4]*(1- clumppped_data[i, 4]))/(2*(clumppped_data[i, 5])^2*clumppped_data[i, 4]*(1- clumppped_data[i, 4]) + ((clumppped_data[i, 6])^2*2*462371*clumppped_data[i, 4]*(1-clumppped_data[i, 4])))
#   f_stats[i] <- (R2[i]*(462371-2))/(1-R2[i])
# }


# # From Katie:
# # 
# # Formula for R^2:
# #f_out[i]<-(2*(beta^2)*(MAF)*(1-(MAF)))/((2*(beta^2)*(MAF)*(1-(MAF)))+((se^2)*(2*N)*(MAF)*(1-(MAF))))

# # -----  Compute F.statistic  -----------------

# # df[i,11] = R2 
# # df[i,9] = N

# #f_out[i]<-(R2(N-2))/(1-R2)

# F.stat<-rep(0,6)
# for(i in 1:6){
#   F.stat[i]<-(df[i,11]*(df[i,9] -2))/(1- df[i,11])
# }

# # Add to table :
# df$F.stat<-F.stat
# #######
##########

##########
# Harmonise and clean up (this merges files as well):
# Current are defaults
action <- as.numeric(args[['--action']])
# action <- 2
print('Harmonising datasets.')
harmonised <- TwoSampleMR::harmonise_data(exposure_dat = exp_data,
                                          outcome_dat = out_data,
                                          action = action # 2=infer positive strand alleles,
                                          # using AF for palindromes
                                          )
# head(harmonised)
# dim(harmonised)

# Save to disk for inspection:
filename <- sprintf('%s.harmonised_tsv',
                    output_file_prefix
                    )
message('Saving harmonised dataset to:')
message(filename)
episcout::epi_write(harmonised, filename)

# # Drop duplicate exposure-outcome summary sets:
# # Drop by sample size, use this if the exposure is binary:
# dat_prune_1 <- power_prune(harmonised, method = 1, dist.outcome = "binary")
# dim(dat_prune_1)
# # Drop by instrument strength:
# dat_prune_2 <- power_prune(harmonised, method = 2, dist.outcome = "binary")
# dim(dat_prune_2)
##########
######################

######################
# Run MR methods

##########
# Run analysis for each SNP individually as default:
pipeMR_single <- function(input = NULL) {
  message('Running MR analysis on each SNP individually.')
  # Obtain MR estimates for each of the selected IVs:
  res_single <- TwoSampleMR::mr_singlesnp(harmonised # Other params:
                                          # parameters = default_parameters(),
                                          # single_method = "mr_wald_ratio",
                                          # all_method = c("mr_ivw",
                                          #                "mr_egger_regression"
                                          #                )
                                          )
  # res_single
  res_single <- TwoSampleMR::generate_odds_ratios(res_single)
  message('Default parameters for TwoSampleMR:')
  message(default_parameters())
  # Save to file:
  filename <- sprintf('%s.results_single_SNP',
                      output_file_prefix
                      )
  message('Saving single SNP analysis to:')
  message(filename)
  episcout::epi_write(res_single, filename)
  # NOTE:
  # empty data.tables get saved as empty files]
  return(res_single)
  }
res_single <- pipeMR_single(input = harmonised)
# res_single
##########

##########
# Run addtional methods:
mr_methods <- as.character(args[['--mr-methods']])
# mr_methods <- 'all'

if (mr_methods == 'none') {
  warning('No additional MR analysis methods requested.')
  } else if (mr_methods == 'main') {
    message('Running MR main methods.')
    mr_main <- TwoSampleMR::mr(harmonised)
    mr_OR_and_CIs <- TwoSampleMR::generate_odds_ratios(mr_main)
    # mr_OR_and_CIs
    } else if (mr_methods == 'all') {
      # Run with more methods:
      message('Running MR all methods.')
      # TwoSampleMR::mr_method_list()
      # TO DO: add single method string option, eg:
      # methods_list <- c("mr_ivw")
      methods_list <- c("mr_wald_ratio", # single SNP only
                        "mr_two_sample_ml",
                        "mr_egger_regression",
                        "mr_egger_regression_bootstrap",
                        "mr_simple_median",
                        "mr_weighted_median",
                        "mr_penalised_weighted_median",
                        "mr_ivw",
                        "mr_ivw_mre",
                        # "mr_ivw_radial", # needs the RadialMR package, run separately
                        "mr_ivw_fe",
                        "mr_simple_mode",
                        "mr_weighted_mode",
                        "mr_weighted_mode_nome",
                        "mr_simple_mode_nome",
                        # "mr_raps", # needs the mr.raps package, run below
                        "mr_sign",
                        "mr_uwr"
                        )
      mr_all <- TwoSampleMR::mr(harmonised, method_list = methods_list)
      # results_mr_all
      # Add OR and CIs:
      mr_OR_and_CIs <- TwoSampleMR::generate_odds_ratios(mr_all)
      # mr_OR_and_CIs
    } else {
      stop('--mr-methods parameter not provided correctly?')
    }

# Save to file:
if (mr_methods == 'main' | mr_methods == 'all') {
  filename <- sprintf('%s.results_all_mr',
                      output_file_prefix
                      )
  message('Saving MR analysis to:')
  message(filename)
  episcout::epi_write(mr_OR_and_CIs, filename)
  }
##########

##########
# Run MR IVW radial:
# args[['--run-radial']] <- TRUE
if (isTRUE(args[['--run-radial']])) {
  message('Running MR IVW radial method.')
  res_radial <- mr(harmonised, method_list = c("mr_ivw_radial"))
  # res_radial
  # Save to file:
  filename <- sprintf('%s.results_mr_radial',
                      output_file_prefix
                      )
  message('Saving MR IVW radial analysis to:')
  message(filename)
  episcout::epi_write(res_radial, filename)
  # # Save full output:
  # needs sink to capture output to stdout from function call though
  # filename <- sprintf('%s.results_mr_radial_full',
  #                     output_file_prefix
  #                     )
  # message('Saving MR IVW radial full output analysis to:')
  # message(filename)
  # cat("MR IVW radial test from RadialMR package\n",
  #     file = filename
  #     )
  # capture.output(res_radial,
  #                file = filename,
  #                append = TRUE
  #                )
  } else {
  message('--run-radial set to FALSE (default)')
    }
##########

##########
# Run MR RAPS:
# args[['--run-raps']] <- TRUE
# TO DO: errors
# mr_raps <- TwoSampleMR::mr(harmonised, method_list = "mr_raps")
if (isTRUE(args[['--run-raps']])) {
  message('Running MR RAPS method.')
  res_raps <- mr(harmonised, method_list = c("mr_raps"))
  # res_raps
  # Save to file:
  filename <- sprintf('%s.results_mr_raps',
                      output_file_prefix
                      )
  message('Saving MR RAPS results to:')
  message(filename)
  episcout::epi_write(res_raps, filename)
  } else {
  message('--run-raps set to FALSE (default)')
    }
##########
######################

######################
# Get more stats

##########
# Convert to MendelianRandomization package to get more tests:
# Largely already provided above
# Errors when passing to MR package
# twoSMR_to_MR_pack <- TwoSampleMR::dat_to_MRInput(dat = harmonised,
#                                                  get_correlations = FALSE
#                                                  )
# Do manually:
twoSMR_to_MR_pack <- MendelianRandomization::mr_input(
  bx = harmonised$beta.exposure,
  bxse = harmonised$se.exposure,
  by = harmonised$beta.outcome,
  byse = harmonised$se.outcome,
  exposure = harmonised$exposure[[1]],
  outcome = harmonised$outcome[[1]],
  snps = harmonised$SNP
  )
# class(twoSMR_to_MR_pack)
# twoSMR_to_MR_pack

# Get I^2 statistic:
# Set up args:
egger_robust <- as.logical(args[['--egger-robust']])
egger_penalized <- as.logical(args[['--egger-penalized']])
egger_distribution <- as.character(args[['--egger-distribution']])
egger_alpha <- as.numeric(args[['--egger-alpha']])
# egger_robust <- FALSE
# egger_penalized <- FALSE
# egger_distribution <- 'normal'
# egger_alpha <- 0.05
# TO DO: check what extra info this provides
egger <- MendelianRandomization::mr_egger(twoSMR_to_MR_pack,
                                          robust = egger_robust, # FALSE
                                          penalized = egger_penalized, # FALSE
                                          distribution = egger_distribution, # "normal"
                                          alpha = egger_alpha # 0.05
                                          )
# class(egger)
# str(egger)

# Save to file:
filename <- sprintf('%s.results_egger_i2',
                    output_file_prefix
                    )
message('Saving Egger analysis to:')
message(filename)
cat("MR Egger test from MendelianRandomization package to obtain I-squared\n",
    file = filename
    )
capture.output(egger,
               file = filename,
               append = TRUE
               )
# NOTE:
# If egger object is NULL, cat will print NULL to file
##########
######################

######################
# Heterogeneity and pleiotropy

##########
# Heterogeneity statistics
hetero_meths <- TwoSampleMR::mr_method_list()
# TO DO: pass as arg:
# hetero_meths <- hetero_meths[which(hetero_meths$heterogeneity_test), 'obj']
hetero_meths <- c("mr_two_sample_ml",
                  "mr_egger_regression",
                  "mr_ivw",
                  # "mr_ivw_radial"#, # requires package ‘RadialMR’
                  "mr_uwr"
                  )

hetero_all <- TwoSampleMR::mr_heterogeneity(harmonised,
                                            method_list = hetero_meths
                                            )
# hetero_all
# TO DO, better calculate as:
#You can use mr_heterogeneity(dat) to obtain Q statistics, and then calculate Isq using
#x <- mr_heterogeneity(dat)
#(x$Q - x$Q_df) / x$Q) * 100

# Save to file:
filename <- sprintf('%s.results_heterogeneity',
                    output_file_prefix
                    )
message('Saving heterogeneity analysis to:')
message(filename)
episcout::epi_write(hetero_all, filename)
# NOTE:
# empty data.tables get saved as empty files
##########

##########
# Horizontal pleiotropy:
pleio <- TwoSampleMR::mr_pleiotropy_test(harmonised)
# pleio
# Save to file:
filename <- sprintf('%s.results_pleiotropy',
                    output_file_prefix
                    )
message('Saving pleiotropy analysis to:')
message(filename)
episcout::epi_write(pleio, filename)
# NOTE:
# empty data.tables get saved as empty files
##########
######################

######################
# Sensitivity analysis

##########
# Obtain MR estimates excluding one IV at a time:
res_loo <- TwoSampleMR::mr_leaveoneout(harmonised)
# parameters = default_parameters(),
# method = mr_ivw
# )
# res_loo
# Save to file:
filename <- sprintf('%s.results_loo',
                    output_file_prefix
                    )
message('Saving leave one out analysis to:')
message(filename)
episcout::epi_write(res_loo, filename)
# NOTE:
# empty data.tables get saved as empty files
##########

##########
# Test that the exposure is upstream of the outcome:
outcome_type <- as.character(args[['--outcome-type']])

run_steiger <- function(outcome_type) {
   out = tryCatch( {
     if ((outcome_type == 'quant') &
         (all(c('pval.exposure',
                'pval.outcome',
                'samplesize.exposure',
                'samplesize.outcome'
                ) %in% names(harmonised)
             )
          )
        ) {
          steiger <- TwoSampleMR::directionality_test(harmonised)
          # Save to file:
          filename <- sprintf('%s.results_steiger',
                              output_file_prefix
                              )
          message('Saving Steiger test to:')
          message(filename)
          episcout::epi_write(steiger, filename)
          } else if (!all(c('pval.exposure',
                            'pval.outcome',
                            'samplesize.exposure',
                            'samplesize.outcome'
                            ) %in% names(harmonised)
                          )
                     ) {
                warning("Missing p-values or sample sizes for exposure and/or outcome, can't calculate.")
              } else if (outcome_type == 'binary') {
                     warning('Steiger test only calculated for quantitative traits in this script.')
  # TO DO:
  # "automated correlations assume quantitative traits. For binary traits
  #  please pre-calculate in r.exposure and r.outcome e.g. using get_r_from_lor()"
  }
},
                  error = function(cond) {
                    message('Directionality test error message: ')
                    message(cond)
                    # Return value in case of error:
                    return(NULL)
                  },
                  warning = function(cond) {
                    message('Directionality test warning message: ')
                    message(cond)
                    # Return value in case of warning:
                    return(NULL)
                  },
                  finally = {
                    message('Directionality test processing done, see results, errors or warnings generated.')
                  }
  )
  return(out)
  }
steiger_results <- run_steiger(outcome_type = outcome_type)
##########

##########
# Run MRPRESSO for outlier detection:
# Input:
MRPRESSO_input <- data.frame(beta_exposure = harmonised$beta.exposure,
                             se_exposure = harmonised$se.exposure,
                             beta_outcome = harmonised$beta.outcome,
                             se_outcome = harmonised$se.outcome,
                             row.names = harmonised$SNP
                             )
# names(MRPRESSO_input)
# MRPRESSO_input

# Set up args:
OUTLIERtest <- as.logical(args[['--OUTLIERtest']]) # TRUE
DISTORTIONtest <- as.logical(args[['--DISTORTIONtest']])# TRUE
NbDistribution <- as.numeric(args[['--NbDistribution']]) # 1000
SignifThreshold <- as.numeric(args[['--SignifThreshold']]) # 0.05

MRPRESSO <- function(OUTLIERtest,
                     DISTORTIONtest,
                     NbDistribution,
                     SignifThreshold
                     ) {
       out = tryCatch( { MRPRESSO::mr_presso(BetaOutcome = "beta_outcome",
                                             BetaExposure = "beta_exposure",
                                             SdOutcome = "se_outcome",
                                             SdExposure = "se_exposure",
                                             OUTLIERtest = OUTLIERtest, # TRUE
                                             DISTORTIONtest = DISTORTIONtest, # TRUE
                                             data = MRPRESSO_input,
                                             NbDistribution = NbDistribution, # 1000
                                             SignifThreshold = SignifThreshold # 0.05
                                             )
                          # MRPRESSO
                          # class(MRPRESSO)
                          # class(MRPRESSO$`Main MR results`)
                          # class(MRPRESSO$`MR-PRESSO results`)
                          # str(MRPRESSO)
                          # TO DO: parse results properly
                          # Save to file:
                          filename <- sprintf('%s.results_presso',
                                              output_file_prefix
                                              )
                          message('Saving MRPRESSO analysis to:')
                          message(filename)
                          cat('MRPRESSO tests from MRPRESSO package to detect pleiotropy (global test) and outliers (outlier test)\n\n',
                          file = filename
                             )
                          capture.output(MRPRESSO,
                                         file = filename,
                                         append = TRUE
                                         )
                          },
                  error = function(cond) {
                    message('MRPRESSO error message: ')
                    message(cond)
                    # Return value in case of error:
                    return(NULL)
                  },
                  warning = function(cond) {
                    message('MRPRESSO warning message: ')
                    message(cond)
                    # Return value in case of warning:
                    return(NULL)
                  },
                  finally = {
                    message('MRPRESSO processing done, see results, errors or warnings generated.')
                  }
  )
  return(out)
  }
mrpresso_results <- MRPRESSO(OUTLIERtest = OUTLIERtest,
                             DISTORTIONtest = DISTORTIONtest,
                             NbDistribution = NbDistribution,
                             SignifThreshold = SignifThreshold
                             )
##########
######################

######################
##########
# TO DO:
# Put results together in a single table
# check if combine_all_mrresults() is helpful
# https://mrcieu.github.io/TwoSampleMR/reference/combine_all_mrresults.html

# # TO DO: will error if any value is empty:
# # add error catch
# all_res <- TwoSampleMR::combine_all_mrresults(#mr_OR_and_CIs,
#                                               res = res_all_OR_and_CIs,
#                                               het = hetero_all,
#                                               plt = pleio,
#                                               sin = res_single,
#                                               ao_slc = F,
#                                               Exp = F,
#                                               split.exposure = F,
#                                               split.outcome = F
#                                               )
#
# # head(all_res[,c("Method","outcome","exposure","nsnp","b","se","pval","intercept","intercept_se","intercept_pval","Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])

# Manually for eg single SNP results:
# echo -e "exposure\toutcome\tid.exposure\tid.outcome\tsamplesize\tSNP\tb\tse\tp\tlo_ci\tup_ci\tor\tor_lci95\tor_uci95" > header_single_SNP.tsv
# cat summary_Wald.tsv | grep -v exposure > summary_Wald.tsv2
# cat header_single_SNP.tsv summary_Wald.tsv2 > summary_Wald.tsv
# rm -f summary_Wald.tsv2
##########
######################

######################
##########
# Additional tests:
# TO DO:
# Unsure what this does:
# TwoSampleMR::enrichment(harmonised) # docs have little info
##########
######################

######################
# Get plots

##########
# Scatter plots:
pipemr_scatter <- function(mr_results = NULL,
                           dat = NULL
                           ) {
  p1 <- TwoSampleMR::mr_scatter_plot(mr_results = mr_OR_and_CIs,
                                     dat = harmonised
                                     )[[1]] + scale_colour_hue()
  # p1
  # class(p1[[1]])
  # class(p1)
  # ggplot2::ggsave('mr_scatter_plot.svg', scale = 1.2)

  # TO DO:
  # add legend eg
  ## Save some text:
  # Methods
  # Legend
  # Interpretation
  # cat(file <- output_file, some_var, '\t', another_var, '\n', append = TRUE)
  message('Generating scatter plot.')
  return(p1)
}

if (mr_methods == 'main' | mr_methods == 'all') {
  p1 <- pipemr_scatter(mr_results = mr_OR_and_CIs,
                       dat = harmonised
                       )
  } else {
    p1 <- NULL
    warning('Scatter plot not generated')
  }
##########

##########
# Get plots
# Single SNP plots:
pipemr_single_SNP_plots <- function(singlesnp_results = NULL,
                                    exponentiate = FALSE,
                                    leaveoneout_results = NULL
                                    ) {
  # Forest plot:
  p2 <- TwoSampleMR::mr_forest_plot(singlesnp_results = res_single,
                                    exponentiate = exponentiate
  )
  # p2
  # ggplot2::ggsave('mr_forest_plot.svg')

  # Leave one out plot:
  p3 <- TwoSampleMR::mr_leaveoneout_plot(leaveoneout_results = res_loo)
  # p3
  # ggplot2::ggsave('mr_leaveoneout_plot.svg')

  # Funnel plot, reciprocal of SE vs MR estimate:
  p4 <- TwoSampleMR::mr_funnel_plot(singlesnp_results = res_single)
  # p4
  # ggplot2::ggsave('mr_funnel_plot.svg')

  # Put together:
  single_SNP_plots <- list(p2[[1]], p3[[1]], p4[[1]])

  # TO DO:
  # add legend
  message('Running MR single SNP plots.')
  return(single_SNP_plots)
}

# Generate plots:
single_plots <- pipemr_single_SNP_plots(singlesnp_results = res_single,
                                        exponentiate = FALSE,
                                        leaveoneout_results = res_loo
                                        )
# single_plots[[3]]

# Save to file:
if (is.null(p1)) {
  grid_plots <- episcout::epi_plots_to_grid(single_plots,
                                            label_size = 12,
                                            align = 'v'
                                            )
} else {
  plot_list <- list(p1, # already extracting plot from TwoSampleMR object
                    single_plots[[1]], # so that it is a plain list
                    single_plots[[2]],
                    single_plots[[3]]
                    )
  grid_plots <- episcout::epi_plots_to_grid(plot_list,
                                            label_size = 12,
                                            align = 'v'
  )
}
  filename <- sprintf('%s.results.svg', # suffix has to be '.svg', '.pdf' etc.
                      output_file_prefix
                      )
  message('Saving plots to:')
  message(filename)
  episcout::epi_plot_cow_save(filename,
                              grid_plots,
                              base_height = 16,
                              base_width = 16
                              )
##########

##########
# TO DO:
# See the 1-to-many forest plot option
# https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html

# and continue with more of the vignette for additional scenarios
##########

##########
# Get an html report:
# mr_report(harmonised)
# This saves a '.md' and '.html' files.
##########
######################

######################
# The end
# Filename to save current R session, data and objects at the end:
if (isTRUE(args[['--session']])) {
	save_session <- sprintf('mr_session_%s.RData', output_file_prefix)
  message(sprintf('Saving an R session image as: %s', save_session))
  save.image(file = save_session, compress = 'gzip')
} else {
  message('Not saving an R session image, this is the default. Specify the --session option otherwise')
}

# If using Rscript and creating plots, Rscript will create the file Rplots.pdf
# by default, it doesn't look like there is an easy way to suppress it, so deleting here:
message('Deleting the file Rplots.pdf...')
system('rm -f Rplots.pdf')
message('Finished successfully.')
sessionInfo()
q()
# Next: run the script for xxx
######################
