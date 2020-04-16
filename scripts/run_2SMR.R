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
TwoSampleMR


Usage and options
=================

Usage: run_2SMR.R (--exposure <exposure_file>) (--outcome <outcome_file>) [options]
       run_2SMR.R [-h | --help]

Options:
  --exposure <exposure_file>      Exposure input file name
  --outcome <outcome_file>        Outcome input file name
  -O <OUTPUT_FILE>                Output file name
  --run-clumping <BOOL>           Clump SNPs with TwoSampleMR [default: TRUE].
  --clump-kb <num>                TwoSampleMR clump_data() parameter [default: 10000].
  --clump-r2 <num>                TwoSampleMR clump_data() parameter [default: 0.001].
  --clump-p1 <num>                TwoSampleMR clump_data() parameter [default: 1].
  --clump-p2 <num>                TwoSampleMR clump_data() parameter [default: 1].
  --action <int>                  TwoSampleMR action to take for harmonising [default: 2].
  --egger-robust <BOOL>           TwoSampleMR mr_egger() parameter [default: FALSE]
  --egger-penalized <BOOL>        TwoSampleMR mr_egger() parameter [default: FALSE]
  --egger-distribution <chr>      TwoSampleMR mr_egger() parameter [default: normal]
  --egger-alpha <num>             TwoSampleMR mr_egger() parameter [default: 0.05]
  --OUTLIERtest <BOOL>            TwoSampleMR mr_egger() parameter [default: TRUE]
  --DISTORTIONtest <BOOL>         TwoSampleMR mr_egger() parameter [default: TRUE]
  --NbDistribution <int>          TwoSampleMR mr_egger() parameter [default: 1000]
  --SignifThreshold <num>         TwoSampleMR mr_egger() parameter [default: 0.05]
  --session                       R session if to be saved [default: NULL]
  -h --help                       Show this screen

Input:

    A tab separated file with headers which conform to the TwoSampleMR
    requirements.

Output:

Various plots and tables for 2SMR.

Requirements:

    library(docopt)
    library(TwoSampleMR)
    library(ggplot2)

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
######################

##########
# TO DO:

# Get proxy SNPs for exposure data first, then match to outcome data, then clump
# (remove weaker instruments in LD) and harmonise
# Currently matching without proxy SNPs


# function to make table from:
# mr_all_BMI_on_CHD.tsv
# single_SNP_wald_BMI_on_CHD.tsv
# heterogeneity_BMI_on_CHD.tsv
# i_squared_BMI_on_CHD.txt
# loo_IVW_BMI_on_CHD.tsv
# mrpresso_BMI_on_CHD.tsv
# pleiotropy_BMI_on_CHD.tsv
# steiger_BMI_on_CHD.tsv

# Other results are:
# mr_all_plots_BMI_on_CHD.svg


# pass as arg:
# TwoSampleMR::mr_method_list() # hardcoded, currently all but radial

# check as Parsa files not consistent in naming
# check clumping Ville did at R2 of 0.1 is correct
##########

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
library(dplyr)
library(ggplot2)
library(cowplot)

library(psych)
library(RadialMR)
library(mr.raps)
library(MRPRESSO)

library(MendelianRandomization)
library(TwoSampleMR)

# source functions from a different R script:
# source(file.path(Rscripts_dir, 'ggtheme.R'))
######################

######################
##########
# Read exposure data:
if (!is.null(args[['--exposure']])) { # for docopt this will be NULL or chr, if boolean
	                            # remove is.null function and test with ==
  exposure_file <- as.character(args[['--exposure']])
  # For tests:
  # setwd('~/xxxx/')
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
  print('You need to provide an exposure input file. This has to be tab separated with headers in the format for TwoSampleMR.')
  stopifnot(!is.null(args[['--exposure']]))
}

print('File being used for exposure data: ')
print(exposure_file)
##########

##########
# Read outcome data:
if (!is.null(args[['--outcome']])) { # for docopt this will be NULL or chr, if boolean
  # remove is.null function and test with ==
  outcome_file <- as.character(args[['--outcome']])
  # For tests:
  # setwd('~/xxxx/')
  # outcome_file <- 'adj_PF_PLTF_WY.2SMR.tsv'

  out_data <- TwoSampleMR::read_outcome_data(filename = outcome_file,
                                             sep = '\t',
                                             phenotype_col = 'Phenotype'
                                             )
  # colnames(out_data)
  # head(out_data)
  # dim(out_data)

  # Stop if arguments not given:
  print('You need to provide an outcome input file. This has to be tab separated with headers in the format for TwoSampleMR.')
  stopifnot(!is.null(args[['--outcome']]))
}

print('File being used for outcome data: ')
print(outcome_file)
##########

##########
# Set output file names:
# suffix <- '2SMR'
if (is.null(args[['-O']])) { # arg is NULL
  # Split infile names at the last '.':
  outfile_name_exp <- strsplit(exposure_file, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  outfile_name_out <- strsplit(outcome_file, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  output_file_prefix <- sprintf('%s_on_%s', outfile_name_exp, outfile_name_out)
  output_file_prefix
  print('Output file prefix not given. Using: ')
  print(output_file_prefix)
} else {
  output_file_prefix <- as.character(args[['-O']])
  # output_file_prefix <- 'testing'
  print(sprintf('Output file prefix provided: %s', output_file_prefix))
  print(output_file_prefix)
}
##########
######################

######################
##########
# First add warning(), stop(), if/else according to number of SNPs
# then move each section as a separate function to pipeMR if needed
# Clump SNPs:
# Currently defaults are set
# Set up args:
clump_kb <- as.numeric(args[['--clump-kb']])
clump_r2 <- as.numeric(args[['--clump-r2']])
clump_p1 <- as.numeric(args[['--clump-p1']])
clump_p2 <- as.numeric(args[['--clump-p2']])

if (isTRUE(args[['--run-clump']])) { # clumping required
  exp_data <- TwoSampleMR::clump_data(dat = exp_data,
                                    clump_kb = clump_kb,
                                    clump_r2 = clump_r2,
                                    clump_p1 = clump_p1,
                                    clump_p2 = clump_p2
                                    )
} else{
  warning('Clumping set to false.')
}
head(exp_data)
dim(exp_data)



# Harmonise and clean up (this merges files as well):
# Current are defaults
action <- as.numeric(args[['--action']])
harmonised <- TwoSampleMR::harmonise_data(exposure_dat = exp_data,
                                          outcome_dat = out_data,
                                          action = action # 2=infer positive strand alleles,
                                          # using AF for palindromes
                                          )
head(harmonised)
dim(harmonised)

# # Drop duplicate exposure-outcome summary sets:
# # Drop by sample size:
# dat_prune_1 <- power_prune(harmonised, method = 1, dist.outcome = "binary")
# dim(dat_prune_1)
# # Drop by instrument strength:
# dat_prune_2 <- power_prune(harmonised, method = 2, dist.outcome = "binary")
# dim(dat_prune_2)
##########

##########
# Run MR and further tests:
results_mr <- TwoSampleMR::mr(harmonised)
results_mr
# default_parameters():
# $test_dist
# [1] "z"
#
# $nboot
# [1] 1000
#
# $Cov
# [1] 0
#
# $penk
# [1] 20
#
# $phi
# [1] 1
#
# $alpha
# [1] 0.05
#
# $Qthresh
# [1] 0.05
#
# $over.dispersion
# [1] TRUE
#
# $loss.function
# [1] "huber"

# Re-run with more methods:
# TwoSampleMR::mr_method_list()
methods_list <- c("mr_wald_ratio", # single SNP only
                  "mr_two_sample_ml",
                  "mr_egger_regression",
                  "mr_egger_regression_bootstrap",
                  "mr_simple_median",
                  "mr_weighted_median",
                  "mr_penalised_weighted_median",
                  "mr_ivw",
                  # "mr_ivw_radial", # needs the RadialMR package, errors
                  "mr_ivw_mre",
                  "mr_ivw_fe",
                  "mr_simple_mode",
                  "mr_weighted_mode",
                  "mr_weighted_mode_nome",
                  "mr_simple_mode_nome",
                  "mr_raps", # needs the mr.raps package
                  "mr_sign",
                  "mr_uwr")

results_mr_all <- TwoSampleMR::mr(harmonised, method_list = methods_list)
results_mr_all
# Add OR and CIs:
res_all_OR_and_CIs <- TwoSampleMR::generate_odds_ratios(results_mr_all)
res_all_OR_and_CIs

# Errors:
# res_radial <- mr(harmonised, method_list = c("mr_ivw_radial"))
# res_radial

# Save to file:
filename <- sprintf('mr_all_%s.txt',
                    output_file_prefix
                    )
print('Saving main MR analysis to:')
print(filename)
episcout::epi_write(res_all_OR_and_CIs, filename)
##########

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
class(twoSMR_to_MR_pack)
twoSMR_to_MR_pack

# Get I^2 statistic:
# Set up args:
egger_robust <- as.character(args[['--egger-robust']])
egger_penalized <- as.character(args[['--egger-penalized']])
egger_distribution <- as.character(args[['--egger-distribution']])
egger_alpha <- as.numeric(args[['--egger-alpha']])

egger <- MendelianRandomization::mr_egger(twoSMR_to_MR_pack,
                                          robust = egger_robust, # FALSE
                                          penalized = egger_penalized, # FALSE
                                          distribution = egger_distribution, # "normal"
                                          alpha = egger_alpha # 0.05
                                          )
class(egger)
str(egger)

# Save to file:
filename <- sprintf('egger_i2_%s.txt',
                    output_file_prefix
                    )
print('Saving Egger analysis to:')
print(filename)
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


##########
# Sensitivity analysis
# Heterogeneity statistics
hetero_meths <- TwoSampleMR::mr_method_list()
hetero_meths <- hetero_meths[which(hetero_meths$heterogeneity_test), 'obj']
hetero_meths <- c("mr_two_sample_ml",
                  "mr_egger_regression",
                  "mr_ivw",
                  # TO DO:
                  # "mr_ivw_radial", # requires package ‘RadialMR’
                  "mr_uwr"
                  )

hetero_all <- TwoSampleMR::mr_heterogeneity(harmonised,
                                            method_list = hetero_meths
                                            )
hetero_all
# Save to file:
filename <- sprintf('heterogeneity_%s.txt',
                    output_file_prefix
                    )
print('Saving heterogeneity analysis to:')
print(filename)
episcout::epi_write(hetero_all, filename)
# NOTE:
# empty data.tables get saved as empty files


# Horizontal pleiotropy:
pleio <- TwoSampleMR::mr_pleiotropy_test(harmonised)
pleio
# Save to file:
filename <- sprintf('pleiotropy_%s.txt',
                    output_file_prefix
                    )
print('Saving pleiotropy analysis to:')
print(filename)
episcout::epi_write(pleio, filename)
# NOTE:
# empty data.tables get saved as empty files


# Obtain MR estimates for each of the selected IVs:
res_single <- TwoSampleMR::mr_singlesnp(harmonised,
                                        parameters = default_parameters(),
                                        single_method = "mr_wald_ratio",
                                        all_method = c("mr_ivw",
                                                       "mr_egger_regression"
                                                       )
                                        )
res_single
# Save to file:
filename <- sprintf('single_SNP_wald_%s.txt',
                    output_file_prefix
                    )
print('Saving single SNP analysis to:')
print(filename)
episcout::epi_write(res_single, filename)
# NOTE:
# empty data.tables get saved as empty files


# Obtain MR estimates excluding one IV at a time:
res_loo <- TwoSampleMR::mr_leaveoneout(harmonised)
# parameters = default_parameters(),
# method = mr_ivw
# )
res_loo
# Save to file:
filename <- sprintf('loo_IVW_%s.txt',
                    output_file_prefix
                    )
print('Saving leave one out analysis to:')
print(filename)
episcout::epi_write(res_loo, filename)
# NOTE:
# empty data.tables get saved as empty files

# Additional tests:
# TO DO:
# Unsure what this does:
TwoSampleMR::enrichment(harmonised)

# Test that the exposure is upstream of the outcome:
if (dim(harmonised)[1] > 1) {
  steiger <- TwoSampleMR::directionality_test(harmonised)
  steiger
  # Save to file:
  filename <- sprintf('steiger_%s.txt',
                      output_file_prefix
                      )
  print('Saving Steiger test to:')
  print(filename)
  episcout::epi_write(steiger, filename)
  } else{
    warning('Too few instruments to run Steiger test, skipping.')
}
##########

##########
# Run MRPRESSO for outlier detection:
# Input:
if (dim(harmonised)[1] > 1) {
  MRPRESSO_input <- data.frame(beta_exposure = harmonised$beta.exposure,
                               se_exposure = harmonised$se.exposure,
                               beta_outcome = harmonised$beta.outcome,
                               se_outcome = harmonised$se.outcome,
                               row.names = harmonised$SNP
                               )
  names(MRPRESSO_input)
  MRPRESSO_input

  # Set up args:
  OUTLIERtest <- as.character(args[['--OUTLIERtest']]) # TRUE
  DISTORTIONtest <- as.character(args[['--DISTORTIONtest']])# TRUE
  NbDistribution <- as.numeric(args[['--NbDistribution']]) # 1000
  SignifThreshold <- as.numeric(args[['--SignifThreshold']]) # 0.05

  MRPRESSO <- MRPRESSO::mr_presso(BetaOutcome = "beta_outcome",
                                  BetaExposure = "beta_exposure",
                                  SdOutcome = "se_outcome",
                                  SdExposure = "se_exposure",
                                  OUTLIERtest = OUTLIERtest, # TRUE
                                  DISTORTIONtest = DISTORTIONtest, # TRUE
                                  data = MRPRESSO_input,
                                  NbDistribution = NbDistribution, # 1000
                                  SignifThreshold = SignifThreshold # 0.05
                                  )
  MRPRESSO
  class(MRPRESSO)
  class(MRPRESSO$`Main MR results`)
  class(MRPRESSO$`MR-PRESSO results`)
  str(MRPRESSO)
  # Save to file:
  filename <- sprintf('mrpresso_%s.txt',
                      output_file_prefix
                      )
  print('Saving MRPRESSO analysis to:')
  print(filename)
  cat('MRPRESSO tests from MRPRESSO package to detect pleiotropy (global test) and outliers (outlier test)\n\n',
      file = filename
    )
  capture.output(MRPRESSO,
                 file = filename,
                 append = TRUE
                 )
  } else{
    warning('Too few instruments to run MRPRESSO test, skipping.')
}
##########

##########
# TO DO:
# Put results together in a single table
# check if combine_all_mrresults() is helpful
# https://mrcieu.github.io/TwoSampleMR/reference/combine_all_mrresults.html

# Output objects from above:
# TwoSampleMR:
# results_mr_all
# res_all_OR_and_CIs
# hetero_all
# pleio
# res_single
# res_loo
# steiger

# MendelianRandomization:
# egger

# Others:
# MRPRESSO

# # TO DO: will error if any value is empty:
# # add error catch
# all_res <- TwoSampleMR::combine_all_mrresults(#results_mr_all,
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
##########

##########
# Get plots
# Scatter plots:
print('Running MR plots.')
p1 <- TwoSampleMR::mr_scatter_plot(results_mr_all, harmonised)[[1]] + scale_colour_hue()
# p1
class(p1[[1]])
class(p1)
# ggplot2::ggsave('mr_scatter_plot.svg', scale = 1.2)

# Forest plot:
p2 <- TwoSampleMR::mr_forest_plot(res_single)
# p2
# ggplot2::ggsave('mr_forest_plot.svg')

# Leave one out plot:
p3 <- TwoSampleMR::mr_leaveoneout_plot(res_loo)
# p3
# ggplot2::ggsave('mr_leaveoneout_plot.svg')

# Funnel plot, reciprocal of SE vs MR estimate:
p4 <- TwoSampleMR::mr_funnel_plot(res_single)
# p4
# ggplot2::ggsave('mr_funnel_plot.svg')

# Put together:
my_plots <- list(p1, p2[[1]], p3[[1]], p4[[1]]) # p1 was already subset
sapply(my_plots, class)
grid_plots <- episcout::epi_plots_to_grid(my_plots,
                                          label_size = 12,
                                          align = 'v'
                                          )
# TO DO:
# add legend
filename <- sprintf('mr_all_plots_%s.svg',
                    output_file_prefix
                    )
print('Saving plots to:')
print(filename)
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
## Save some text:
# Methods
# Legend
# Interpretation
# cat(file <- output_file, some_var, '\t', another_var, '\n', append = TRUE)
######################

######################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(xxx))
#objects_to_save <- (c('xxx_var'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# Filename to save current R session, data and objects at the end:
if (!is.null(args[['--session']])) { # arg is NULL
	save_session <- sprintf('%s_%s.RData', output_file_prefix)
  print(sprintf('Saving an R session image as: %s', save_session))
  save.image(file = save_session, compress = 'gzip')
} else {
  print('Not saving an R session image, this is the default. Specify the --session option otherwise')
}

# If using Rscript and creating plots, Rscript will create the file Rplots.pdf
# by default, it doesn't look like there is an easy way to suppress it, so deleting here:
print('Deleting the file Rplots.pdf...')
system('rm -f Rplots.pdf')
print('Finished successfully.')
sessionInfo()
q()

# Next: run the script for xxx
######################
