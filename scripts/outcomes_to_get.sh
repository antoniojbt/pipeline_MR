#!/usr/bin/env bash

# run individual bash commands to extract data for each pmid

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
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O SLE_26502338 --pmid 26502338
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O MS_31604244 --pmid 31604244
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O T1D_25751624 --pmid 25751624
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O RA_24390342 --pmid 24390342
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O psoriasis_28537254 --pmid 28537254
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O shared_paeds_AID_26301688 --pmid 26301688
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O CD_Trynka_2011_22057235 --pmid 22057235
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O Ellinghaus_shared_5_AID_2017_26974007 --pmid 26974007
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O IBD_28067908 --pmid 28067908
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O AS_23749187 --pmid 23749187
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O AS_20062062 --pmid 20062062
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O AS_21743469 --pmid 21743469
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O IBD_shared_26192919 --pmid 26192919
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O UC_21297633 --pmid 21297633
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O UC_20228799 --pmid 20228799
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O vitiligo_27723757 --pmid 27723757
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O vitiligo_22561518 --pmid 22561518
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O thyroid_22922229 --pmid 22922229
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O Behcet_27548383 --pmid 27548383
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O PAD_31285632 --pmid 31285632
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O atherosclerosis_21909108 --pmid 21909108
./get_ieugwasr.R --instruments cis_pqtl_instruments_SNPs_only.txt -O atherosclerosis_30510157 --pmid 30510157
###########################
