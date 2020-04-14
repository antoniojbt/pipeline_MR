##########
# Cytokines and blood cell traits MR
# March 2020
# Antonio
# Run MR for Ville cytokines (exposures) with Parsa et al blood cells traits (outcomes)
##########

##########
# TO DO:

# TO DO: check
# Get proxy SNPs for exposure data first, then match to outcome data, then clump
# (remove weaker instruments in LD) and harmonise
# Currently matching without proxy SNPs

# TO DO:
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
##########

##########
# setwd('~/Documents/quickstart_projects/projects/cytokines_blood_traits_MR/results/results_MR/')
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
library(dplyr)
library(ggplot2)
library(cowplot)

library(psych)
library(RadialMR)
library(mr.raps)
library(MRPRESSO)

library(MendelianRandomization)
library(TwoSampleMR)
##########

##########
# Read in files
# TO DO: pass as arg:
# exposure_file <- 'cis_pqtl_instruments.2SMR.tsv'
exposure_file <- as.character(args[1])

# outcome_file <- 'pqtl_in_adj_PF_PLTF_WY.2SMR.tsv'
# outcome_file <- 'adj_PF_PLTF_WY.2SMR.tsv'
outcome_file <- as.character(args[2])

# Get exposure instruments and match to outcome data:
exp_data <- TwoSampleMR::read_exposure_data(filename = exposure_file,
                                            sep = '\t',
                                            clump = F,
                                            phenotype_col = 'Phenotype'
                                            )
colnames(exp_data)
head(exp_data)
dim(exp_data)

out_data <- TwoSampleMR::read_outcome_data(filename = outcome_file,
                                           sep = '\t',
                                           phenotype_col = 'Phenotype'
                                           )
colnames(out_data)
head(out_data)
dim(out_data)
##########

##########
# Set output file names:
# TO DO: check as Parsa files not consistent in naming
outfile_name_exp <- strsplit(exposure_file, split = '[.]')[[1]][1]
outfile_name_exp <- strsplit(outfile_name_exp, split = '_')
outfile_name_exp <- paste(outfile_name_exp[[1]][1],
                          outfile_name_exp[[1]][2],
                          sep = '_'
                          )
outfile_name_exp
outfile_name_out <- strsplit(outcome_file, split = '[.]')[[1]][1]
outfile_name_out

# TO DO: pass as arg:
# Others depending on how much flexibility to give:
# TwoSampleMR::clump_data() # more important
# TwoSampleMR::harmonise_data() 'action' arg # more important
# TwoSampleMR::mr_method_list() # hardcoded, currently all but radial
# MendelianRandomization::mr_egger() # save defaults to legend or log
# MRPRESSO::mr_presso # save defaults to legend or log
##########

##########
# Clump SNPs:
# TO DO: pass as args:
# Currently defaults are set
# TO DO: check clumping is needed
exp_data <- TwoSampleMR::clump_data(dat = exp_data,
                                    clump_kb = 10000,
                                    clump_r2 = 0.001,
                                    clump_p1 = 1,
                                    clump_p2 = 1
                                    )
head(exp_data)
dim(exp_data)


# Harmonise and clean up (this merges files as well):
# TO DO: pass as args:
# Current are defaults
harmonised <- TwoSampleMR::harmonise_data(exposure_dat = exp_data,
                                   outcome_dat = out_data,
                                   action = 2 # 2=infer positive strand alleles,
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
                  "mr_uwr"
                  )

results_mr_all <- TwoSampleMR::mr(harmonised, method_list = methods_list)
results_mr_all
# Add OR and CIs:
res_all_OR_and_CIs <- TwoSampleMR::generate_odds_ratios(results_mr_all)
res_all_OR_and_CIs

# Errors:
# res_radial <- mr(harmonised, method_list = c("mr_ivw_radial"))
# res_radial

# Save to file:
filename <- sprintf('mr_all_%s_on_%s.txt', outfile_name_exp, outfile_name_out)
filename
episcout::epi_write(res_all_OR_and_CIs, filename)
##########

##########
# Convert to MendelianRandomization package to get more tests:
# Largely already provided above
# TO DO: get args
# TO DO: errors when passing to MR package
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
egger <- MendelianRandomization::mr_egger(twoSMR_to_MR_pack,
                                          robust = FALSE,
                                          penalized = FALSE,
                                          distribution = "normal",
                                          alpha = 0.05
                                          )
class(egger)
str(egger)

# Save to file:
filename <- sprintf('egger_i_squared_%s_on_%s.txt',
                    outfile_name_exp,
                    outfile_name_out
                    )
filename
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
filename <- sprintf('heterogeneity_%s_on_%s.txt',
                    outfile_name_exp,
                    outfile_name_out
                    )
filename
episcout::epi_write(hetero_all, filename)
# NOTE:
# empty data.tables get saved as empty files


# Horizontal pleiotropy:
pleio <- TwoSampleMR::mr_pleiotropy_test(harmonised)
pleio
# Save to file:
filename <- sprintf('pleiotropy_%s_on_%s.txt',
                    outfile_name_exp,
                    outfile_name_out
                    )
filename
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
filename <- sprintf('single_SNP_wald_%s_on_%s.txt',
                    outfile_name_exp,
                    outfile_name_out
                    )
filename
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
filename <- sprintf('loo_IVW_%s_on_%s.txt',
                    outfile_name_exp,
                    outfile_name_out
                    )
filename
episcout::epi_write(res_loo, filename)
# NOTE:
# empty data.tables get saved as empty files

# Additional tests:
# TO DO:
# Unsure what this does:
TwoSampleMR::enrichment(harmonised)

# TO DO: add error catch
# # Test that the exposure is upstream of the outcome:
# # TO DO: will error if not enough instruments
# steiger <- TwoSampleMR::directionality_test(harmonised)
# steiger
# # Save to file:
# filename <- sprintf('steiger_%s_on_%s.txt',
#                     outfile_name_exp,
#                     outfile_name_out
#                     )
# filename
# episcout::epi_write(steiger, filename)
# # NOTE:
# # empty data.tables get saved as empty files
##########

##########
# TO DO: add error catch
# # Run MRPRESSO for outlier detection:
# # Input:
# MRPRESSO_input <- data.frame(beta_exposure = harmonised$beta.exposure,
#                              se_exposure = harmonised$se.exposure,
#                              beta_outcome = harmonised$beta.outcome,
#                              se_outcome = harmonised$se.outcome,
#                              row.names = harmonised$SNP
#                              )
# names(MRPRESSO_input)
# MRPRESSO_input
#
# # TO DO: will error if not enough instruments
# MRPRESSO <- MRPRESSO::mr_presso(BetaOutcome = "beta_outcome",
#                                 BetaExposure = "beta_exposure",
#                                 SdOutcome = "se_outcome",
#                                 SdExposure = "se_exposure",
#                                 OUTLIERtest = TRUE,
#                                 DISTORTIONtest = TRUE,
#                                 data = MRPRESSO_input,
#                                 NbDistribution = 1000,
#                                 SignifThreshold = 0.05
#                                 )
# MRPRESSO
# class(MRPRESSO)
# class(MRPRESSO$`Main MR results`)
# class(MRPRESSO$`MR-PRESSO results`)
# str(MRPRESSO)
# # Save to file:
# filename <- sprintf('mrpresso_%s_on_%s.txt',
#                     outfile_name_exp,
#                     outfile_name_out
#                     )
# filename
# cat('MRPRESSO tests from MRPRESSO package to detect pleiotropy (global test) and outliers (outlier test)\n\n',
#     file = filename
#     )
# capture.output(MRPRESSO,
#                file = filename,
#                append = TRUE
#                )
##########

##########
# TO DO:
# Put results together in a single table
# TO DO:
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

# TO DO: add error catch
# # TO DO: will error if any value is empty:
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
filename <- sprintf('mr_all_plots_%s_on_%s.svg',
                    outfile_name_exp,
                    outfile_name_out
                    )
filename
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

##########
sessionInfo()
q()
##########
