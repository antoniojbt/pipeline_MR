#!/usr/bin/env bash

# Needs a file called headers.txt (output from get_headers.sh)
# Process headers and add a phenotype column so that files match TwoSampleMR input.

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
PHENOTYPE=$2

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

# Add Phenotype column, remove .tsv.gz from column info:
use_GNU_cmds ${FILE} | sed -r "s/$/\t${PHENOTYPE}/" | sed '1d' | cat <(echo -e $(tail -n 1 headers.txt)) - | gzip > ${FILE}.2SMR_tsv.gz &&
    echo 'Finished processing:' ${FILE}