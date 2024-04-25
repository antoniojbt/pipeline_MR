#!/usr/bin/env bash

# Get headers from GWAS summary files
# Output needs manual processing, including changing the delimiter for '\t' Then pass to header_to_2SMR.sh
# so that echo -e can pass the string correctly
# In vim this would need e.g. :%s/^I/\\t/gc
# where ^I is the tab insert and needs to be done with ctrl-v then 
# key so that the literal binding is inserted
# Delimeters and headers vary by file so easier if done manually

# TwoSampleMR read_exposure_data() requirements are:
#https://mrcieu.github.io/TwoSampleMR/reference/read_exposure_data.html
# phenotype_col = "Phenotype"
# snp_col = "SNP"
# beta_col = "beta"
# se_col = "se"
# eaf_col = "eaf"
# effect_allele_col = "effect_allele"
# other_allele_col = "other_allele"
# pval_col = "pval"
# units_col = "units"
# ncase_col = "ncase"
# ncontrol_col = "ncontrol"
# samplesize_col = "samplesize"
# gene_col = "gene"
# id_col = "id"

###########################
# Set bash script options

# exit when a command fails:
set -o errexit

# exit if any pipe commands fail:
set -o pipefail

# exit if there are undeclared variables:
set -o nounset

# trace what gets executed:
set -o xtrace
set -o errtrace
###########################

###########################

FILE=$1

# Avoid Mac's zcat:
cmd_name=gzcat
use_GNU_cmds() {
    if ! cmd_loc="$(type -p "$cmd_name")" || [[ -z $cmd_loc ]]; then
    #if type "gzcat" 2>/dev/null; then
        #alias zcat=gzcat
        ${cmd_name} "$@"
    else
        ${cmd_name} "$@"
    fi
}

# Check conda env being used:
echo $CONDA_PREFIX

# Get two lines of the header to keep a record of the original one, modify the second header to match
# TwoSampleMR input requirements
cat <(head -n 1 <(use_GNU_cmds $FILE)) <(head -n 1 <(use_GNU_cmds $FILE)) > headers.txt
