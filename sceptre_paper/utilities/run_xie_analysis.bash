#!/bin/bash

# This script runs the Xie data analysis. It is assumed that this script is being executed from within the utilities directory.

# Set the machine and number of processors.
machine=local
n_processors=20
precomputation_complete=FALSE

echo Running on machine $machine with $n_processors processors.

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

echo Check the availability of the required packages
Rscript $code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"check_packages_0.R $code_dir

# Make sure the sceptre and katsevich2020 packages are up-to-date.
# bash build_and_install_package.bash sceptre $machine
# bash build_and_install_package.bash katsevich2020 $machine

echo Initialize the offsite directory structure.
# Rscript $code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"check_directory_structure_1.R $code_dir $offsite_dir

echo Download the data. Note that one of the downloads must be done manually.
# Rscript $code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"download_data_2.R $code_dir $offsite_dir

echo Pre-process the data.
# Rscript $code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"pre_process_data_3.R $code_dir $offsite_dir

echo Construct model covariate matrix and perform quality control.
# Rscript $code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"quality_control_4.R $code_dir $offsite_dir

echo Compute gene-gRNA pairs to analyze.
# Rscript $code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"select_gRNA_gene_pair_4.1.R $code_dir $offsite_dir

# Locate the parameter file
parameter_file=$code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"sceptre_function_args.R

echo Run the sceptre analysis at scale.
sceptre_at_scale_bash_dir=$code_dir"/sceptre_package/sceptre_at_scale"
bash $sceptre_at_scale_bash_dir/"run_sceptre_at_scale.bash" $sceptre_at_scale_bash_dir $offsite_dir $parameter_file $n_processors $precomputation_complete

echo Run the negative binomial regression analysis at scale.
nb_at_scale_bash_dir=$code_dir"/sceptre_paper/"nb_regression_at_scale
bash $nb_at_scale_bash_dir/"run_nb_regression_at_scale.bash" $nb_at_scale_bash_dir $offsite_dir $parameter_file $n_processors
wait

echo Preprocess the sceptre results for downstream analysis.
# Rscript $code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"append_simple_names_results_5.R $code_dir $offsite_dir

echo Run bulk RNA-seq analysis.
# Rscript $code_dir"/sceptre_paper/analysis_drivers/analysis_drivers_xie/"bulk_validation_6.R $code_dir $offsite_dir
