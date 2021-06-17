#!/bin/bash

# This script runs the entire analysis, starting with initializing the directory structure and downloading the required packages and ending with creating the figures. To run this script on a new computer, update the "local" file paths in the utilities/get_file_paths.bash file. Change "/Users/timbarry/Box/SCEPTRE-manuscript/SCEPTRE" to the location of the cloned github repository on the new computer, and change "/Volumes/tims_new_drive/research/sceptre_files" to an empty directory on the new computer to store all the results (preferably on an external hard drive).
# Note 1: It is assumed that this script is being executed from within the sceptre_paper/utilities directory.
# Note 2: This script likely would take several months to run on a standard laptop. We ran the majority of the analysis on a high-performance computer cluster.

# Run analysis on Gasperini data
bash run_gasperini_analysis.bash

# Run analysis on Xie data
bash run_xie_analysis.bash

# Run simulation study
bash run_simulation.bash

# Create the figures
bash make_figures.bash
