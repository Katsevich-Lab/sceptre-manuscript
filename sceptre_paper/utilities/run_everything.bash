#!/bin/bash

# This script runs the entire analysis (from creating the directory structure, checking package availability, and downloading the raw data to running the analysis and making plots) from start to finish. To run this script on a new computer, update the "local" file paths in the get_file_paths.bash file. To do this, change "/Users/timbarry/Box/SCEPTRE/SCEPTRE" to the location of the cloned github repository on the new computer, and change "/Volumes/tims_new_drive/research/sceptre_files" to an empty directory on the new computer (preferably on an external hard drive).
# Note: One of the datasets for this analysis must be downloaded manually. Check sceptre_paper/analysis_drivers/analysis_drivers_xie/download_data_2.R for details.
# Note: It is assumed that this script is being executed from within the sceptre_paper/utilities directory.

# Run analysis on Gasperini data
bash run_gasperini_analysis.bash

# Run analysis on Xie data
bash run_xie_analysis.bash

# Run simulation study
bash run_simulation.bash

# Create the figures
bash make_figures.bash