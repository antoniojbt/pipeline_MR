#!/usr/bin/env python3

##############
# Get all the modules needed
# System:
import os
import subprocess

# Directory:
#os.chdir('/Users/antoniob/Documents/quickstart_projects/data/external/MR_data/parsa_flow_cytometry')

def create_symlinks(file_to_read, path, var1, var2):
    '''
    Create symlinks for Parsa files to keep
    '''
    # 
    with open(file_to_read) as f:
        next(f) # skip header
        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            #print(fields)
            full_path = os.path.join(path, fields[var1])
            #print(full_path)
            #print(var1, var2)
            proc = subprocess.run(["ln", "-s", full_path, "."])
            # If renaming symlinks:
            # new_name = var2
            #proc = subprocess.run(["ln", "-s", full_path, var2])
            #print("Command is {}".format(proc.args))
    return()

create_symlinks(file_to_read = 'parsa_traits_to_keep.tsv',
                path = '/Users/antoniob/Documents/quickstart_projects/data/external/MR_data/parsa_flow_cytometry/parsa_cytokines_summary_stats_2_28_apr_2020/',
                var1 = 2, # should be "filename"
                var2 = 0 # should be "Standard Abbreviation"
                )
