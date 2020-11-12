#!/bin/bash

# This bash script runs the Gasperini data analysis. The code is commented to increase ease of adoption and use. It is assumed that this bash file is being executed from within the utilities directory.

# Set the machine.
machine=local
n_processors=20

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

echo Check the availability of the required packages
Rscript $code_dir"/analysis_drivers_gasp/"check_packages_0.R $code_dir

# Make sure the sceptre and katsevich2020 packages are up-to-date.
bash build_and_install_package.bash sceptre $machine
bash build_and_install_package.bash katsevich2020 $machine

echo Initialize the offsite directory structure.
# Rscript $code_dir"/analysis_drivers_gasp/"check_directory_structure_1.R $code_dir $offsite_dir

echo Download the data.
# Rscript $code_dir"/analysis_drivers_gasp/"download_raw_data_2.R $code_dir $offsite_dir

echo Pre-process the data.
# Rscript $code_dir"/analysis_drivers_gasp/"pre_process_data_3.R $code_dir $offsite_dir

echo Perform quality control.
# Rscript $code_dir"/analysis_drivers_gasp/"perform_quality_control_4.R $code_dir $offsite_dir

# Locate the parameter file
parameter_file=$code_dir"/analysis_drivers_gasp/sceptre_function_args.R"

echo Run the sceptre analysis at scale.
sceptre_at_scale_bash_dir=$code_dir"/functions_at_scale/sceptre_at_scale"
bash $sceptre_at_scale_bash_dir/"run_sceptre_at_scale.bash" $sceptre_at_scale_bash_dir $offsite_dir $parameter_file $n_processors

echo Run the negative binomial regression at scale.
nb_at_scale_bash_dir=$code_dir"/functions_at_scale/nb_regression_at_scale"
bash $nb_at_scale_bash_dir/"run_nb_regression_at_scale.bash" $nb_at_scale_bash_dir $offsite_dir $parameter_file $n_processors
