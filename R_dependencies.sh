#!/usr/bin/env bash

# Install packages for pipeline_MR and run_2SMR.R
# Requires a conda environment setup with tidyverse, python 3.0 and R 4


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
conda install -y docopt \
                 r-docopt \
                 r-data.table \
                 r-cowplot

conda install -y r-svglite

R -e "install.packages('devtools', repos = 'https://cloud.r-project.org')"
R -e "install.packages('psych', repos = 'https://cloud.r-project.org')"
R -e "install.packages('MendelianRandomization', repos = 'https://cloud.r-project.org')"

R -e "library(devtools) ; install_github('AntonioJBT/episcout')"
R -e "library(devtools) ; install_github('MRCIEU/TwoSampleMR')"
R -e "library(devtools) ; install_github('mrcieu/ieugwasr')"

R -e "library(devtools) ; install_github('qingyuanzhao/mr.raps')"

# These get installed by packages above:
#R -e "library(devtools) ; install_github('rondolab/MR-PRESSO')"
#R -e "library(devtools) ; install_github('WSpiller/RadialMR')"


echo 'Done installing'

#R -e "install.packages('', repos = 'https://cloud.r-project.org')"
#R -e "library(devtools) ; install_github('')"
###########################