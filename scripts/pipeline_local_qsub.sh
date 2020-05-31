#!/bin/sh

#####
# PBS directives need to go first

#Job parameters, throughput:
#PBS -lwalltime=02:00:00
#PBS -lselect=1:ncpus=32:mem=62gb

# Copy current terminal environment to nodes
# If using conda activate first in submission terminal
#PBS -V
# this can be checked with qstat -f job_ID and checking the 
# Variable_List description in the output
# Some environment variables are copied across regardless though
# such as PATH and LANG

# Email alert at (b)eginning, (e)nd and (a)bortion of execution
# Probably not available though
##PBS -m bea
##PBS -M aberlang@ic.ac.uk
#####

#####
# Treat as normal bash script, set bash options:
# exit when a command fails:
set -o errexit
# exit if any pipe commands fail:
set -o pipefail
# exit if there are undeclared variables:
set -o nounset
# trace what gets executed:
set -o xtrace
set -o errtrace
#####

#####
# Copy any input files needed across to nodes:
cd ${PBS_O_WORKDIR}
#####

#####
# Actual commands of interest to run:
conda activate ds
python pipeline_MR.py make full -v 5 --local
#####

#####
# Copy results back to working directory from nodes:
#cp ${OUTFILE} $PBS_O_WORKDIR
#####
