#!/usr/bin/env bash

# run individual bash commands for 2SMR for each pmid extracted as outcome

###########################
# Set bash script options

# exit when a command fails:
#set -o errexit

# exit if any pipe commands fail:
set -o pipefail

# exit if there are undeclared variables:
set -o nounset

# trace what gets executed:
set -o xtrace
set -o errtrace
###########################

###########################
# Run commands:
../run_2SMR.R --exposure ../cis_pqtl_instruments.2SMR.tsv --stop-clumping \
              --outcome ../outcome_data_available/AS_23749187.tsv \
              -O pqtl_AS_23749187


../run_2SMR.R --exposure ../cis_pqtl_instruments.2SMR.tsv --stop-clumping \
              --outcome ../outcome_data_available/IBD_28067908.tsv \
              -O pqtl_IBD_28067908


../run_2SMR.R --exposure ../cis_pqtl_instruments.2SMR.tsv --stop-clumping \
              --outcome ../outcome_data_available/IBD_shared_26192919.tsv \
              -O pqtl_IBD_shared_26192919


../run_2SMR.R --exposure ../cis_pqtl_instruments.2SMR.tsv --stop-clumping \
              --outcome ../outcome_data_available/RA_24390342.tsv \
              -O pqtl_RA_24390342


../run_2SMR.R --exposure ../cis_pqtl_instruments.2SMR.tsv --stop-clumping \
              --outcome ../outcome_data_available/SLE_26502338.tsv \
              -O pqtl_SLE_26502338


../run_2SMR.R --exposure ../cis_pqtl_instruments.2SMR.tsv --stop-clumping \
              --outcome ../outcome_data_available/T1D_25751624.tsv \
              -O pqtl_T1D_25751624


../run_2SMR.R --exposure ../cis_pqtl_instruments.2SMR.tsv --stop-clumping \
              --outcome ../outcome_data_available/UC_21297633.tsv \
              -O pqtl_UC_21297633
