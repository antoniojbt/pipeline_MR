'''
pipeline_name
=============

:Author: |author_name|
:Release: |version|
:Date: |today|


Overview
========

|long_description|


Purpose
=======

.. briefly describe the main purpose and methods of the pipeline



Usage and options
=================

These are based on CGATPipelines_ and Ruffus_, not docopt.

.. _CGATPipelines: https://github.com/CGATOxford/CGATPipelines

.. _Ruffus: http://www.ruffus.org.uk/


For command line help type:

    python pipeline_pq_example.py --help

Configuration
=============

This pipeline is built using a Ruffus/CGAT approach. You need to have Python,
Ruffus, CGAT core tools and any other specific dependencies needed for this
script.

A configuration file was created at the same time as this script.

Use this to extract any arbitrary parameters that could be changed in future
re-runs of the pipeline.


Input files
===========

.. Describe the input files needed, urls for reference and preferably place
example data somewhere.


Pipeline output
===============

.. Describe output files and results


Requirements
============

cgat-core as well as the following software need to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* R >= 3.4
* Python >= 3.5

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
    ''' Get the absolute path to where this function resides. Useful for
    determining the user's path to a package. If a sub-directory is given it
    will be added to the path returned. Use '..' to go up directory levels. '''
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


@transform('*.csv', # adj_PF_PLTF_WZ.csv
           suffix('.csv'),
           '.csv_to_tsv'
           )
def csv_to_tsv(infile, outfile):
    '''
    Convert simple (no quotes) comma separated files to tab separated
    '''
    # for parsa csv files:
    statement = '''cat %(infile)s | tr ',' '\\t' > %(outfile)s'''

    P.run(statement)


#@follows(csv_to_tsv)
@transform('*.csv_to_tsv',
           suffix('.csv_to_tsv'),
           '.rg_csv_to_tsv',
           'exposure_instruments.txt', # cis_pqtl_instruments_SNPs_only.txt
           )
def grep_SNPs(infile1, outfile, infile2):
    '''
    Grep exposure SNPs from outcome files using ripgrep
    '''
    # rg -wf cis_pqtl_instruments_SNPs_only.txt \
    # adj_PF_PLTF_WZ.tsv > \
    # pqtl_in_adj_PF_PLTF_WZ.tsv
    statement = '''rg -wf %(infile2)s %(infile1)s > %(outfile)s '''
    
    P.run(statement)


#@follows(grep_SNPs)
@transform('*.rg_csv_to_tsv',
           suffix('.rg_csv_to_tsv'),
           '.out_2SMR_tsv'
           )
def parsa_to_2SMR(infile, outfile):
    '''
    Change column names to TwoSampleMR format for Parsa et al files 
    '''
    
    statement = '''Rscript parsa_to_2SMR.R %(infile)s %(outfile)s
                '''

    P.run(statement)


# Only needs to be run once and is specific to one file:
#@follows(grep_SNPs)
@transform('cis_*.csv', # cis_pqtl_instruments.csv
           suffix('.csv'),
           '.exp_2SMR_tsv' # cis_pqtl.exp_2SMR_tsv
           )
def ville_to_2SMR(infile, outfile):
    '''
    Change column names to TwoSampleMR format for Ville pQTL file
    '''
    
    statement = '''Rscript ville_to_2SMR.R %(infile)s %(outfile)s
                '''

    P.run(statement)


#@follows(parsa_to_2SMR)
@transform('*.out_2SMR_tsv',
           suffix('.out_2SMR_tsv'),
           '.out_2SMR_touch', # multiple outputs
           'cis_pqtl.exp_2SMR_tsv' # use a generic name
           )
def ville_parsa_MR(outcome, outfile, exposure):
    '''
    Run MR analysis for Ville cytokines (exposure) and Parsa et al blood
    cell traits (outcome)
    '''
    
    statement = '''Rscript ville_parsa_MR.R %(exposure)s %(outcome)s ; 
                   touch %(outfile)s
                '''

    P.run(statement)


#@follows(ville_parsa_MR)
@transform('*.svg',
           suffix('.svg'),
           '.pdf'
           )
def svg_to_pdf(infile, outfile):
    '''
    Convert svg plots to PDFs and merge into single document
    '''
    #ln -s $(grep -L 'Insufficient' ../results_all/%(infile)s) . ;
    statement = ''' 
                    rsvg-convert -a -f pdf -o %(outfile)s %(infile)s
                '''

    P.run(statement)


#@follows(svg_to_pdf)
@merge(svg_to_pdf, 'all_plots.pdf')
def combine_pdfs(infiles, summary_file):
    '''
    Combine individual PDFs into single file
    '''
#pdfunite %(infiles)s %(summary_file)s
    statement = ''' pdfunite *.pdf %(summary_file)s
                '''

    P.run(statement)


# TO DO: combine results

################

################
# Copy to log enviroment from conda:
@follows(ville_parsa_MR)
@originate('conda_info.txt')
def conda_info(outfile):
    '''
    Save to logs conda information and packages installed.
    '''
    packages = 'conda_packages.txt'
    channels = 'conda_channels.txt'
    environment = 'environment.yml'

    statement = '''conda info -a > %(outfile)s ;
                   conda list -e > %(packages)s ;
                   conda list --show-channel-urls > %(channels)s ;
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
    ''' Generates html and pdf versions of restructuredText files
        using sphinx-quickstart pre-configured files (conf.py and Makefile).
        Pre-configured files need to be in a pre-existing report directory.
        Existing reports are overwritten.
    '''
    report_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               'pipeline_report'
                                               ))
    print('Copying report templates from: {}'.format(report_path))

    if (os.path.exists(report_dir) and
            os.path.isdir(report_dir) and not
            os.listdir(report_dir)):
        statement = '''cp %(report_path)s/* pipeline_report ;
                       cd {} ;
                       ln -s ../pipeline.yml . ;
                       make html ;
                       ln -sf _build/html/report_pipeline_pq_example.html . ;
                       make latexpdf ;
                       ln -sf _build/latex/pq_example.pdf .
                    '''.format(report_dir)
        E.info('''Building pdf and html versions of your rst files in
                  {}.'''.format(report_dir))
        P.run(statement)

    elif (os.path.exists(report_dir) and
            os.path.isdir(report_dir) and
            os.listdir(report_dir)):
        sys.exit(''' {1} exists, not overwriting. You can manually run:
                       cd {1} ;
                       ln -s ../pipeline.yml . ;
                       make html ;
                       ln -sf _build/html/report_XXXX.html . ;
                       make latexpdf ;
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
#                    ln -s pipeline_report/_build/html/index.hmtl %(project_name)s.html ;
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
