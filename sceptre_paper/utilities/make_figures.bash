#!/bin/bash

# This script creates the main components of Figures 1-4, S1-S3.

# Set the machine.
machine=uberduo

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

# Check packages and directory structure.
Rscript $code_dir"/sceptre_paper/plotting/check_packages_and_dirs_0.R" $code_dir

# Create the plots
figures=( 1 2 3 4 5 S1 S2 S3 S4)
for i in "${figures[@]}"
do
  Rscript $code_dir"/sceptre_paper/plotting/Figure"$i".R" $code_dir $offsite_dir
done
