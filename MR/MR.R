#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
script_name
===============

Author: |author_names|
Release: |version|
Date: |today|


Purpose
=======

|description|


Usage and options
=================

These are based on docopt_ for R:

https://github.com/docopt/docopt.R
https://cran.r-project.org/web/packages/docopt/index.html

To run, type:
    Rscript script_name -I <INPUT_FILE> [options]

Usage: script_name (-I <INPUT_FILE>)
       script_name [options]
       script_name [-h | --help]

Options:
  -I <INPUT_FILE>                 Input file name
  -O <OUTPUT_FILE>                Output file name
  --session                       R session if to be saved
  -h --help                       Show this screen
  -var                            some numeric argument [default: 0.001].

Input:

    A tab separated file with headers. This is read with data.table and stringsAsFactors = FALSE

Output:



Requirements:

    library(docopt)
    library(data.table)
    library(ggplot2)

Documentation
=============

    For more information see:

    |url|
' -> doc

# Load docopt:
library(docopt, quietly = TRUE)
# Retrieve the command-line arguments:
args <- docopt(doc)
# See:
# https://cran.r-project.org/web/packages/docopt/docopt.pdf
# docopt(doc, args = commandArgs(TRUE), name = NULL, help = TRUE,
# version = NULL, strict = FALSE, strip_names = !strict,
# quoted_args = !strict)

# Print to screen:
str(args)
# Within the script specify options as:
# args[['--session']]
# args $ `-I` == TRUE
######################

######################
# Logging
# This can be taken care of by CGAT Experiment.py if running as a pipeline.
# Otherwise there seem to be few good alternatives. A workaround is the code in:
# XXXX/project_quickstart/templates/script_templates/logging.R
# It does not run on its own though, needs copy/pasting for now.
######################

######################
# Load a previous R session, data and objects:
#load('R_session_saved_image_order_and_match.RData', verbose=T)
######################

######################
# This function allows other R scripts to obtain the path to a script directory
# (ie where this script lives). Useful when using source('some_script.R')
# without having to pre-specify the location of where that script is.
# This is taken directly from:
# How to source another_file.R from within your R script molgenis/molgenis-pipelines Wiki
# https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
# Couldn't find a licence at the time (12 June 2018)
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }

    if (!is.null(this.file)) return(dirname(this.file))

    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from = 1, length.out = length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))

    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
Rscripts_dir <- LocationOfThisScript()
print('Location where this script lives:')
Rscripts_dir
# R scripts sourced with source() have to be in the same directory as this one
# (or the path constructed appropriately with file.path) eg:
#source(file.path(Rscripts_dir, 'moveme.R')) #, chdir = TRUE)
######################

######################
# Import libraries
# source('http://bioconductor.org/biocLite.R')
# biocLite()
library(ggplot2)
library(data.table)
library(svglite) # prefer over base R svg()
# https://github.com/Rdatatable/data.table/wiki/Getting-started
# source functions from a different R script:
#source(file.path(Rscripts_dir, 'moveme.R')) #, chdir = TRUE)
######################

######################
##########
# Read files, this is with data.table:
if (!is.null(args[['-I']])) { # for docopt this will be NULL or chr, if boolean
	                            # remove is.null function and test with ==
  input_name <- as.character(args[['-I']])
  # For tests:
  # input_name <- 'XXX'
  # setwd('~/xxxx/')
  input_data <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers.')
  stopifnot(!is.null(args[['-I']]))
}

print('File being used: ')
print(input_name)
##########

##########
# Set output file names:
suffix <- 'my_output'
if (is.null(args[['-O']])) { # arg is NULL
  # Split infile name at the last '.':
  output_name <- strsplit(input_name, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  output_file_name <- sprintf('%s.%s', input_name, suffix)
  print('Output file name not given. Using: ')
  print(output_file_name)
} else {
  output_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  output_file_name <- sprintf('%s.%s', output_name, suffix)
  print(sprintf('Output file name provided: %s', output_file_name))
  print(output_file_name)
}
##########
######################

######################
# Explore data:
class(input_data)
dim(input_data) # nrow(), ncol()
str(input_data)
input_data # equivalent to head() and tail()
setkey(input_data) # memory efficient and fast
key(input_data)
tables()
colnames(input_data)
input_data[, 2, with = FALSE] # by column position, preferable by column name to avoid silent bugs
######################

######################
# Process variables
str(input_data)
# The following is much simpler with a dataframe (R base) instead of a data.table
# What you use depends on memory, preference, etc.
# Convert a data.table to a data.frame:
# setDF(input_data) # changes by reference
input_data_df <- as.data.frame(input_data) # create a new object
class(input_data_df)
class(input_data)
# Pass variables as parameters for data.table
# If re-using, automating, etc. you might want to do this,
# so you can then specify variables fom the command line.
# Generally it's just easier (and safer and more readable) to name variables explicitely
# for project specific analysis.
# https://stackoverflow.com/questions/10675182/in-r-data-table-how-do-i-pass-variable-parameters-to-an-expression?rq=1
# Create a function:
pass_var_dt <- function(var_name) {
  a_var <- as.character(var_name)
  a_var_p <- parse(text = a_var)
}
x_var_name <- 'XXX' # still needed
x_var <- pass_var_dt('XXX')
head(input_data[, eval(x_var)])
# Convert variables in data.table
# https://stackoverflow.com/questions/7813578/convert-column-classes-in-data-table
x_var_factor <- as.factor(input_data[, eval(x_var)])
input_data[, c(x_var_name):= x_var_factor, with = F]
str(input_data)
str(input_data_df)

# Setup more variables:
y_var_name <- 'YYY'
y_var <- pass_var_dt('YYY')
head(input_data[, eval(y_var)])
var3_name <- 'age'
var3 <- pass_var_dt('age')
head(input_data[, eval(var3)])
######################

######################
# What's the question?
# What's the hypothesis?
# Descriptive:
nrow(input_data)
length(which(complete.cases(input_data) == TRUE))
summary(input_data)
summary(input_data[, c(2:3)])

# Get the mean for one var specified above:
input_data[, .(mean = mean(eval(var3), na.rm = TRUE))] # drop with and put column name, usually better practice
# Get the mean for all columns except the first one in data.table:
input_data[, lapply(.SD, mean), .SDcols = c(2:ncol(input_data))]
# Specify columns to get some summary stats (numeric variables):
colnames(input_data)
cols_summary <- c(3, 5, 6)
# Control precision for printing, can be nightmarish. Here enforce printing e.g. 0.00
# See:
# https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r
# Modified here so that FUN is function to run, x the number to format
# and k the number of decimals to show.
# This function could/should be moved to a separate script and sourced here.
specify_decimal <- function(FUN, x, k) trimws(format(round(FUN(x, na.rm = TRUE), k), nsmall = k))

# data.table with summary stats:
desc_stats <- input_data[, sapply(.SD, function(x) c(mean = specify_decimal(mean, x, 2),
                                                     median = specify_decimal(median, x, 2),
                                                     SD = specify_decimal(sd, x, 2),
                                                     min = specify_decimal(min, x, 2),
                                                     max = specify_decimal(max, x, 2)
                                                     # quantile_25 = quantile(x, 0.25, na.rm = TRUE), # will error for non-numeric
                                                     # quantile_75 = quantile(x, 0.75, na.rm = TRUE) # will error for non-numeric
                                                     )),
                         .SDcols = cols_summary]
desc_stats

# Make a table from this:
desc_stats <- as.data.frame(desc_stats)
class(desc_stats)
# Add rownames as column with label:
desc_stats$statistics <- rownames(desc_stats)
# Re-order columns:
colnames(desc_stats)
desc_stats <- desc_stats[, c("statistics", "age", "glucose", "BMI")]
desc_stats

# Save file:
fwrite(desc_stats, output_file_name,
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
######################

######################
# Some basic exploratory plots
# Plot here or use a separate plot_template.R script (preferable if processing
# large datasets, process first, save, plot separately):
plot_name <- svg(sprintf('%s_%s_%s.svg', input_name, x_var_name, y_var_name))
# par(mfrow = c(1, 3)) # rows, cols
boxplot(input_data_df[,  y_var_name] ~ input_data_df[,  x_var_name])
dev.off()
######################

######################
####
# Histogram overlaid with kernel density curve
# http://www.cookbook-r.com/Graphs/Plotting_distributions_(ggplot2)/

# Setup:
plot_name <- sprintf('%s_%s_histogram.svg', input_name, x_var)
x_var_label <- x_var
# Plot:
ggplot(input_data, aes(x = input_data[, x_var])) +
       geom_histogram(aes( y = ..density..), # Histogram with density instead of count on y-axis
                 binwidth = 0.5,
                 colour = "black", fill = "white") +
       geom_density(alpha = 0.2, fill = "#FF6666") + # Overlay with transparent density plot
       ylab('density') +
       xlab(x_var_label)
# Save to file:
ggsave(plot_name)
# Prevent Rplots.pdf from being generated. ggsave() without weight/height opens a device.
# Rscript also saves Rplots.pdf by default, these are deleted at the end of this script.
dev.off()
####

####
# A boxplot
# Setup:
var3_label <- var3
var3_factor <- factor(input_data[, var3])
y_var_label <- y_var
plot_name <- sprintf('%s_%s_%s_boxplot.svg', input_name, var3, y_var)
# Plot:
ggplot(input_data, aes(x = var3_factor, y = input_data[, y_var], fill = var3_factor)) +
       geom_boxplot() +
       ylab(y_var_label) +
       xlab(var3_label) +
       theme_classic() +
       theme(legend.position = 'none')
# Save to file:
ggsave(plot_name)
# Prevent Rplots.pdf from being generated. ggsave() without weight/height opens a device.
# Rscript also saves Rplots.pdf by default, these are deleted at the end of this script.
dev.off()
####

####
# Scatterplot and legend:
# http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
# Setup:
x_var_label <- x_var
y_var_label <- y_var
plot_name <- sprintf('%s_%s_%s_scatterplot.svg', input_name, x_var, y_var)
# Plot:
ggplot(input_data, aes(x = input_data[, x_var], y = input_data[, y_var], colour = var3_factor)) +
       geom_point() +
       geom_smooth(method = lm) +
       ylab(y_var_label) +
       xlab(x_var_label) +
       labs(colour = var3) +
       theme_classic()
# Save:
ggsave(plot_name)
# Prevent Rplots.pdf from being generated. ggsave() without weight/height opens a device.
# Rscript also saves Rplots.pdf by default, these are deleted at the end of this script.
dev.off()
####
######################

######################
# # Run in parallel blurb
# # Define cores:
# num_cores <- max(1, detectCores() - 1)
# # Create output holder:
# list_length <- length(unique(df[[some_col]]))
# # List of unique individuals to extract:
# list_unique_ids <- unique(df[[some_id]])
# # Output holder:
# unique_inds <- vector(mode = 'list', length = list_length)
# # Loop through in parallel:
# doFuture::registerDoFuture()
# future::availableCores()
# future::availableWorkers()
# future::plan(multiprocess, workers = num_cores)
# my_out <- foreach(i = 1:length(list_unique_ids)) %dopar% {
# 	# Do something here:
# 	#
# }
# gc() # Clear memory, may not be necessary, especially if running with functions
######################

######################
# Some inferential stats
# Make a table from the output of the linear regression
######################

######################
## Save some text:
# Methods
# Legend
# Interpretation
# cat(file <- output_file, some_var, '\t', another_var, '\n', append = TRUE)
######################

######################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(xxx))
#objects_to_save <- (c('xxx_var'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# Filename to save current R session, data and objects at the end:
if (!is.null(args[['--session']])) { # arg is NULL
	save_session <- sprintf('%s_%s.RData', output_name, suffix)
  print(sprintf('Saving an R session image as: %s', save_session))
  save.image(file = save_session, compress = 'gzip')
} else {
  print('Not saving an R session image, this is the default. Specify the --session option otherwise')
}

# If using Rscript and creating plots, Rscript will create the file Rplots.pdf
# by default, it doesn't look like there is an easy way to suppress it, so deleting here:
print('Deleting the file Rplots.pdf...')
system('rm -f Rplots.pdf')
print('Finished successfully.')
sessionInfo()
q()

# Next: run the script for xxx
######################
