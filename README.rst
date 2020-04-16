.. include:: substitution_vars.rst

.. GitHub doe not render rst substitutions

.. copy across your travis "build..." logo so that it appears in your Github page

.. .. image:: https://travis-ci.org/|github_user|/|project_name|.svg?branch=master
    :target: https://travis-ci.org/|github_user|/|project_name|

.. do the same for ReadtheDocs image:

.. note that if your project is called project_Super readthedocs will convert
.. it to project-super

.. .. image:: https://readthedocs.org/projects/|project_name|/badge/?version=latest
    :target: http://|project_name|.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

 .. Edit manually:

.. .. Zenodo gives a number instead, this needs to be put in manually here:
   .. image:: https://zenodo.org/badge/#######.svg
      :target: https://zenodo.org/badge/latestdoi/#####

**IN PROGRESS**


################################################
pipeline_MR
################################################


.. The following is a modified template from RTD
    http://www.writethedocs.org/guide/writing/beginners-guide-to-docs/#id1

.. For a discussion/approach see 
    http://tom.preston-werner.com/2010/08/23/readme-driven-development.html

cgatcore/Ruffus pipeline for Mendelian randomisation analysis based mainly on the R package TwoSampleMR.

Allows easier running and reporting when needing to process many exposures and/or outcomes.

Input needed
-------------

- Minimum input in a tsv file for both exposure and outcome data using TwoSample MR headers:
  + SNP (rsID)
  + Effect Allele
  + Beta (effect size)
  + SE

- Additional columns (get externally if not provided):
  + Chromosome
  + Position
  + Effect allele frequency (preferably from study sample)
  + Other Allele

- Additional columns (manually added):
  + p-value (from GWAS)
  + Sample size (from GWAS)
  + Exposure phenotype name for plot labels
  + Outcome phenotype name for plot labels


Outputs
--------

Various MR plots and tables

Requirements
------------

Main requirements:

* cgatcore
* R >= 3.6
* Python >= 3.6
* and various R packages including:
	episcout
	dplyr
	ggplot2
	cowplot
	psych
	RadialMR
	mr.raps
	MRPRESSO
	MendelianRandomization
	TwoSampleMR


Installation
------------

.. code-block:: bash
   
    pip install git+git://github.com/AntonioJBT/pipeline_MR.git

Dependencies need to be installed manually for now.


To use
------

.. code-block:: bash

    # Download test files, e.g.:
    wget -nH -np -r --cut-dirs=4 -A .tsv https://github.com/AntonioJBT/pipeline_MR/tree/master/tests/
    python pipeline_MR --help
    python pipeline_MR printconfig
    python pipeline_MR show full
    python pipeline_MR make full -v 2 -p 4 --local
     

Contribute
----------

- `Issue Tracker`_
  
.. _`Issue Tracker`: github.com/AntonioJBT/pipeline_MR/issues

- `Source Code`_
  
.. _`Source Code`: github.com/AntonioJBT/pipeline_MR

- Pull requests welcome!


Support
-------

If you have any issues, pull requests, etc. please report them in the issue tracker. 


