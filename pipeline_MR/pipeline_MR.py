'''
pipeline_MR
=============

:Author: |author_name|
:Release: |version|
:Date: |today|


Overview
========

|long_description|


Purpose
=======

Run two sample Mendelian Randomisation analysis.



Usage and options
=================

These are based on cgatcore and ruffus.

.. _cgatcore: https://github.com/cgat-developers/cgat-core

.. _ruffus: https://github.com/cgat-developers/ruffus


For command line help type:

    pipeline_MR --help
    pipeline_MR show full
    pipeline_MR printconfig
    pipeline_MR config # creates a copy of the yaml file in the working directory
    pipeline_MR make full --local # run locally
    nohup pipeline_MR make full & # submit to cluster, needs DRMAA
    # Check the outputs in output_summary

Configuration
=============

This pipeline is built using a Ruffus/CGAT approach. You need to have Python, Ruffus, CGAT core tools and any other specific dependencies needed for this script.

A configuration file was created at the same time as this script.

Use this to extract any arbitrary parameters that could be changed in future re-runs of the pipeline.


Input files
===========

Requires exposure and outcome tab-separated files as required by the R package TwoSampleMR. Additionally:
  - Outcome file names must be gzipped and have the suffix ".out_2SMR_tsv.gz"
  - The exposure file name (can contain multiple exposures) must be called "exposure.2SMR_tsv" and be decompressed.
  - Only one exposure file is accepted but many outcome files can be passed
  - To have phentoype names in plots and tables include a 'Phenotype' column in the exposure file
  - The outcome and exposure files must at least contain:
    + SNP (rsID)
    + Effect Allele
    + Other Allele
    + Beta (effect size)
    + SE
- A single column file with no header called "exposure_instruments.txt" containing all exposure SNPs must be present in the working directory


Pipeline output
===============

Plots and tables with 2SMR


Requirements
============

cgat-core as well as the following software need to be in the path:

* R
* Various TwoSampleMR R packages

Documentation
=============

    For more information see:

        |url|

'''
################
# Get modules needed:
import sys
import os
import re
import subprocess
import glob

# Pipeline:
from ruffus import *

# Database:
import sqlite3

# CGAT tools:
import cgatcore.iotools as iotools
import cgatcore.pipeline as P
import cgatcore.experiment as E

# For gzipped files:
import binascii
import gzip
################

################
# Get locations of source code (this file)
    # os.path.join note: a subsequent argument with an '/' discards anything
    # before it
    # For function to search path see: 
    # http://stackoverflow.com/questions/4519127/setuptools-package-data-folder-location
# MANIFEST.in file instructs the project_quickstart/templates folder to be included in installation

_ROOT = os.path.abspath(os.path.dirname(__file__))
def getDir(path = _ROOT):
    ''' Get the absolute path to where this function resides. Useful for determining the user's path to a package. If a sub-directory is given it will be added to the path returned. Use '..' to go up directory levels.
    '''
   # src_top_dir = os.path.abspath(os.path.join(_ROOT, '..'))
    src_dir = _ROOT
    return(os.path.abspath(os.path.join(src_dir, path)))
################

################
# Load options from the config file
P.get_parameters(
        ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
            "../pipeline.yml",
            "pipeline.yml"],
        )

PARAMS = P.PARAMS
################

################
# Specific pipeline tasks
# Tools called need the full path or be directly callable

#print('Checking ripgrep')
#import subprocess
#print(subprocess.run(['which', 'rg']))


# TO DO: sort out task completion checking

#####
# Helper functions:
# Unzip files if needed:
def gunzip_files(infile):
    '''
    Check if files are compressed ('.gz') and uncompress.
    '''

    with open(infile, 'rb') as f:
        # gzip files start with this:
        if binascii.hexlify(f.read(2)) == b'1f8b':
            f = gzip.open(f, 'rb')
            file_content = f.read()
            f.close()
        else:
            E.info('File is not gzip compressed')
    return(file_content)


# Zip files if needed:
def gzip_files(infile):
    '''
    Check if files are uncompressed ('.gz') and compress.
    '''
    with open(infile, 'rb') as f:
        # gzip files start with this:
        if binascii.hexlify(f.read(2)) == b'1f8b':
            E.info('File is already gzip compressed')
        else:
            # compress a gzip file
            input = open(infile, 'rb')
            s = input.read()
            input.close()
            
            outfile = '{}.gz'.format(infile)
            output = gzip.GzipFile(outfile, 'wb')
            output.write(s)
            output.close()
    return


# Remove strings from ruffus infiles when more than one is passed as it will be verbatim
def remove_str_ruffus_files(infiles):
    ''' Remove strings which get passed verbatim in ruffus '''
    infiles = str(infiles)
    infiles = infiles.replace("'", "")
    infiles = infiles.replace("[", "")
    infiles = infiles.replace("]", "")
    infiles = infiles.replace(",", "")
    
    return(infiles)


# Keep as utility function, call manually though:
@transform('*.csv',
           suffix('.csv'),
           '.csv_to_tsv'
           )
def csv_to_tsv(infile, outfile):
    '''
    Convert simple (no quotes) comma separated files to tab separated.
    '''
    # Gunzip if compressed:
    #infile = gunzip_files(infile)

    # for parsa csv files:
    statement = '''cat %(infile)s | tr ',' '\\t' > %(outfile)s'''

    P.run(statement)
#####


#####
@transform('*.out_2SMR_tsv.gz', # exposure and outcome files must have this suffix to  be picked up
           suffix('.out_2SMR_tsv.gz'),
           '.rg_2SMR_tsv',
           'exposure_instruments.txt',
           )
def grep_SNPs(infile1, outfile, infile2):
    '''
    Grep exposure SNPs from outcome files using ripgrep.
    '''
    # Gunzip if compressed:
    #infile = gunzip_files(infile1)
    
    # statement = ''' gunzip -c %(infile1)s | head -n 1 > %(infile1)s.header &&
    #                 rg -wf %(infile2)s <(gunzip -c %(infile1)s) | cat %(infile1)s.header - > %(outfile)s && rm -f %(infile1)s.header
    #             '''

    statement = ''' head -n 1 <(zcat < %(infile1)s) > %(infile1)s.header &&
                    rg -wf %(infile2)s <(zcat < %(infile1)s) | cat %(infile1)s.header - > %(outfile)s && rm -f %(infile1)s.header
                '''

    P.run(statement)


@follows(grep_SNPs)
@transform('*.rg_2SMR_tsv',
           suffix('.rg_2SMR_tsv'),
           '.2SMR_touch', # multiple outputs
           'exposure.2SMR_tsv') # must have this name
def run_2SMR(outcome, outfile, exposure):
    '''
    Run MR analysis.
    '''

    run_2SMR_options = PARAMS['run_2SMR']
    if run_2SMR_options == None:
        run_2SMR_options = ''
    else:
        pass

    statement = '''Rscript run_2SMR.R --exposure %(exposure)s \
                                      --outcome %(outcome)s \
                                      %(run_2SMR_options)s &&
                   touch %(outfile)s
                '''

    P.run(statement)


@follows(mkdir('output_summary'), run_2SMR)
@merge('*.results.svg',
       'svgs_with_results.txt'
       )
def get_svgs(infiles, outfile):
    '''
    Get svg plots that have results.
    '''

    # Convert infiles list to str and remove characters:
    infiles = remove_str_ruffus_files(infiles)

    # || : is the same as "or pass"; place the cmd in parenthesis before piping
    # 'true' will also return exit code 0 but can be a builtin in some shells
    # https://unix.stackexchange.com/questions/330660/prevent-grep-from-exiting-in-case-of-nomatch
    statement = ''' (grep -L "Insufficient" %(infiles)s || :) > %(outfile)s &&
                    cat %(outfile)s | while read line ; do ln -rs ${line} -t output_summary/ ; done
                '''

    P.run(statement)


@follows(get_svgs)
@transform('output_summary/*.results.svg',
           suffix('.results.svg'),
           '.results.pdf'
           )
def svg_to_pdf(infile, outfile):
    '''
    Convert svg plots to PDFs and merge into a single document.
    '''
    jobs_limit = 1 # looks like sometimes PDFs get damaged with multiple
                   # processors running and writing in the same dir
    statement = ''' rsvg-convert -a -f pdf -o %(outfile)s %(infile)s
                '''

    P.run(statement)


@follows(svg_to_pdf)
@merge(svg_to_pdf, 'output_summary/all_plots.pdf')
def combine_pdfs(infiles, summary_file):
    '''
    Combine individual PDFs into a single file.
    ''' 

    if (len(infiles) > 1):
        # Convert infiles list to str and remove characters:
        infiles = remove_str_ruffus_files(infiles)
        
        statement = ''' pdfunite %(infiles)s %(summary_file)s
                    '''

        P.run(statement)                

    else:
        E.warn('Only one pdf file present, not combining files.')


@follows(combine_pdfs)
@merge('*.results_single_SNP',
       'mr_single_SNP_summary.tsv'
       )
def single_SNP_summary(infiles, outfile):
    '''
    Combine all results from single SNP analysis into a single table.
    '''
    # Convert infiles list to str and remove characters:
    infiles = remove_str_ruffus_files(infiles)
    
    statement = ''' cat %(infiles)s | grep -v -e exposure -e Inverse -e Egger > %(outfile)s &&
                    echo -e "exposure\\toutcome\\tid.exposure\\tid.outcome\\tsamplesize\\tSNP\\tb\\tse\\tp" | cat - %(outfile)s > %(outfile)s2 &&
                    mv -f %(outfile)s2 %(outfile)s &&
                    mv %(outfile)s output_summary/
                '''

    P.run(statement)


@follows(single_SNP_summary)
@merge('*.results_all_mr',
       'mr_summary.tsv'
       )
def main_mr_summary(infiles, outfile):
    '''
    Combine all results from the main MR analyses into a single table.
    '''
    # Convert infiles list to str and remove characters:
    infiles = remove_str_ruffus_files(infiles)

    statement = ''' cat %(infiles)s | grep -v id.exposure > %(outfile)s &&
                    echo -e "id.exposure\\tid.outcome\\toutcome\\texposure\\tmethod\\tnsnp\\tb\\tse\\tpval\\tlo_ci\\tup_ci\\tor\\tor_lci95\\tor_uci95" | cat - %(outfile)s > %(outfile)s2 &&
                    mv -f %(outfile)s2 %(outfile)s &&
                    mv %(outfile)s output_summary/
                '''

    P.run(statement)


@follows(main_mr_summary)
@transform('output_summary/mr_summary.tsv',
           suffix('mr_summary.tsv'),
           'IVW_only.tsv'
          )
def IVW_summary(infile, outfile):
    '''
    Extract IVW results and sort them by p-value.
    '''
    # TO DO: add FDR
    statement = ''' cat <(head -n 1 %(infile)s) <(grep Inverse %(infile)s) | sort -t$'\\t' -k9 -g > %(outfile)s
                '''

    P.run(statement)


# TO DO: collect MR Radial if called for
# TO DO: collect MR RAPS if called for


@follows(IVW_summary)
@merge('*.results_heterogeneity',
       'heterogeneity_summary.tsv'
       )
def heterogeneity_summary(infiles, outfile):
    '''
    Combine all results from the heterogeneity analysis into a single table.
    '''
    # Convert infiles list to str and remove characters:
    infiles = remove_str_ruffus_files(infiles)

    statement = ''' cat %(infiles)s | grep -v id.exposure > %(outfile)s &&
                    echo -e "id.exposure\\tid.outcome\\toutcome\\texposure\\tmethod\\tQ\\tQ_df\\tQ_pval" | cat - %(outfile)s > %(outfile)s2 &&
                    mv -f %(outfile)s2 %(outfile)s &&
                    mv %(outfile)s output_summary/
                '''

    P.run(statement)


# # TO DO:
# # mr_egger_i2_pqtl_xxx.txt
# # This needs parsing to get exposure name from filename and results from free text file, see eg:
# cat ../results_all/egger_i_squared_cis_pqtl_on_DF_EOSI_WX.txt 
# filename from run_2SMR.R is %s.results_egger_i2


@follows(heterogeneity_summary)
@merge('*.results_loo',
       'loo_summary.tsv'
       )
def loo_summary(infiles, outfile):
    '''
    Combine all results from the leave one out analysis into a single table.
    '''
    # Convert infiles list to str and remove characters:
    infiles = remove_str_ruffus_files(infiles)

    statement = ''' cat %(infiles)s | grep -v exposure > %(outfile)s &&
                    echo -e "exposure\\toutcome\\tid.exposure\\tid.outcome\\tsamplesize\\tSNP\\tb\\tse\\tp" | cat - %(outfile)s > %(outfile)s2 &&
                    mv -f %(outfile)s2 %(outfile)s &&
                    mv %(outfile)s output_summary/
                '''

    P.run(statement)

# TO DO: PRESSO results need to be parsed in run_2SMR
@follows(loo_summary)
@merge('*.results_presso',
       'presso_summary.tsv'
       )
def presso_summary(infiles, outfile):
    '''
    Combine all results from the MR PRESSO analysis into a single table.
    '''
    # Convert infiles list to str and remove characters:
    infiles = remove_str_ruffus_files(infiles)

    statement = ''' cat %(infiles)s | grep -v exposure > %(outfile)s &&
                    echo -e "exposure\\toutcome\\tid.exposure\\tid.outcome\\tsamplesize\\tSNP\\tb\\tse\\tp" | cat - %(outfile)s > %(outfile)s2 &&
                    mv -f %(outfile)s2 %(outfile)s &&
                    mv %(outfile)s output_summary/
                '''

    P.run(statement)


@follows(presso_summary)
@merge('*.results_pleiotropy',
       'pleiotropy_summary.tsv'
       )
def pleiotropy_summary(infiles, outfile):
    '''
    Combine all results from the pleiotropy analysis into a single table.
    '''
    # Convert infiles list to str and remove characters:
    infiles = remove_str_ruffus_files(infiles)

    statement = ''' cat %(infiles)s | grep -v exposure > %(outfile)s &&
                    echo -e "id.exposure\\tid.outcome\\toutcome\\texposure\\tegger_intercept\\tse\\tpval" | cat - %(outfile)s > %(outfile)s2 &&
                    mv -f %(outfile)s2 %(outfile)s &&
                    mv %(outfile)s output_summary/
                '''

    P.run(statement)


@follows(pleiotropy_summary)
@merge('*.results_steiger',
       'steiger_summary.tsv'
       )
def steiger_summary(infiles, outfile):
    '''
    Combine all results from the Steiger directional analysis into a single table.
    '''
    # Convert infiles list to str and remove characters:
    infiles = remove_str_ruffus_files(infiles)

    if len(glob.glob('*.results_steiger')) >= 1 :
        statement = ''' cat %(infiles)s | grep -v exposure > %(outfile)s &&
                        echo -e "id.exposure\\tid.outcome\\texposure\\toutcome\\tsnp_r2.exposure\\tsnp_r2.outcome\\tcorrect_causal_direction\\tsteiger_pval" | cat - %(outfile)s > %(outfile)s2 &&
                        mv -f %(outfile)s2 %(outfile)s &&
                        mv %(outfile)s output_summary/
                    '''
        P.run(statement)

    else:
        E.warn('No Steiger results to process.')

    return
#####
################

################
# Clean up files:
@follows(mkdir('2SMR_analysis'),
         mkdir('ripgrepped_SNPs'),
         steiger_summary
         )
def cleanup():
    '''
    Re-organise files to have separate directories
    '''
    
    statement = ''' mv *results* 2SMR_analysis/ &&
                    mv *.rg_2SMR_tsv ripgrepped_SNPs/ &&
                    rm -f output_summary/*.results.svg &&
                    rm -f output_summary/*.results.pdf
                '''

    P.run(statement)
################

################
# Copy to log enviroment from conda:
@follows(cleanup)
@originate('conda_info.txt')
def conda_info(outfile):
    '''
    Save to logs conda information and packages installed.
    '''
    packages = 'conda_packages.txt'
    channels = 'conda_channels.txt'
    environment = 'environment.yml'

    statement = '''conda info -a > %(outfile)s &&
                   conda list -e > %(packages)s &&
                   conda list --show-channel-urls > %(channels)s &&
                   conda env export > %(environment)s
                '''
    P.run(statement)
################

################
# Create the "full" pipeline target to run all functions specified
@follows(conda_info)
@originate('pipeline_complete.touch')
def full(outfile):
    statement = 'touch %(outfile)s'
    P.run(statement)
################

################
# Build report with pre-configured files using sphinx-quickstart
# Convert any svg files to PDF if needed:
@transform('*.svg', suffix('.svg'), '.pdf')
def svg_to_pdf2(infile, outfile):
    '''
    Simple conversion of svg to pdf files with inkscape
    '''
    statement = '''
                inkscape --without-gui \
                         --export-area-drawing \
                         --export-margin=2 \
                         --file=%(infile)s \
                         --export-pdf=%(outfile)s
                '''
    P.run(statement)


# Build the report:
report_dir = 'pipeline_report'
@follows(svg_to_pdf)
@follows(mkdir(report_dir))
def make_report():
    ''' Generates html and pdf versions of restructuredText files using sphinx-quickstart pre-configured files (conf.py and Makefile). Pre-configured files need to be in a pre-existing report directory. Existing reports are overwritten.
    '''
    report_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               'pipeline_report'
                                               ))
    print('Copying report templates from: {}'.format(report_path))

    if (os.path.exists(report_dir) and
            os.path.isdir(report_dir) and not
            os.listdir(report_dir)):
        statement = '''cp %(report_path)s/* pipeline_report &&
                       cd {} &&
                       ln -s ../pipeline.yml . &&
                       make html &&
                       ln -sf _build/html/report_pipeline_pq_example.html . &&
                       make latexpdf &&
                       ln -sf _build/latex/pq_example.pdf .
                    '''.format(report_dir)
        E.info('''Building pdf and html versions of your rst files in
                  {}.'''.format(report_dir))
        P.run(statement)

    elif (os.path.exists(report_dir) and
            os.path.isdir(report_dir) and
            os.listdir(report_dir)):
        sys.exit(''' {1} exists, not overwriting. You can manually run:
                       cd {1} &&
                       ln -s ../pipeline.yml . &&
                       make html &&
                       ln -sf _build/html/report_XXXX.html . &&
                       make latexpdf &&
                       ln -sf _build/latex/XXXX.pdf .
                       Or delete the folder and re-run make_report
                 '''.format(report_dir)
                 )

    else:
        sys.exit(''' The directory "pipeline_report" does not exist.
                     Are the paths correct?
                     Template files were tried to be copied from:
                     {}
                     You can also manually copy files and run "make html" or
                     "make latexpdf".
                 '''.format(report_path)
                 )

    return

#    if (os.path.exists('pipeline_report/_build/html/index.hmtl') and
#       os.path.exists(os.path.join('pipeline_report/_build/latex/',
#                                   project_name, '.pdf'))):
#        statement = '''
#                    ln -s pipeline_report/_build/html/index.hmtl %(project_name)s.html &&
#                    ln -s pipeline_report/_build/latex/%(project_name)s.pdf .
#                    '''
#        E.info('''Done, links to the pdf and html versions of your rst files are in the main
#               folder.''')
#        P.run()
#
#    else:
#        E.info('''
#               The html and/or latex/pdf files did not build correctly. See the
#               logs and go into pipeline_report to find out. You can also try
#               building the report manually with make html and make latexpdf.
#               ''')
#        sys.exit()
################
# This is if pipeline is called directly from CLI:
# TO DO:
# Check if docopt and argparse can play to show my_pipeline options and P.py
# options

def main():
    sys.exit(P.main(sys.argv))

#def main(argv=None):
#    if argv is None:
#        argv = sys.argv
#    P.main(argv)
################

################
# Otherwise just end pipeline as normally (with sys.exit commented):
if __name__ == "__main__":
    main()
    #sys.exit(P.main(sys.argv))
################
