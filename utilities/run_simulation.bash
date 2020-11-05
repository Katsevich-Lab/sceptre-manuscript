#!/bin/bash

# This bash script runs the simulation analysis. It is assumed that the script is being executed from the utilities directory.

# Set the machine.
machine=local

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

# Make sure the sceptre and katsevich2020 packages are up-to-date.
bash build_and_install_package.bash sceptre $machine
bash build_and_install_package.bash katsevich2020 $machine

# The simulation assumes that the Xie data already has been downloaded. Ensure all simulation folders have been created.
Rscript $code_dir"/analysis_drivers_xie/"check_directory_structure_1.R $code_dir $offsite_dir

# Locate the Xie and simulation parameter files.
parameter_file_xie=$code_dir"/analysis_drivers_xie/sceptre_function_args.R"
parameter_file_sim=$code_dir"/simulations/simulation_parameters.R"

echo Extract a high-quality gRNA and gene on which to base the simulation.
Rscript $code_dir"/simulations/obtain_data_1.R" $offsite_dir $parameter_file_xie

echo Generate the simulation data.
Rscript $code_dir"/simulations/fit_model_generate_data_2.R" $offsite_dir $parameter_file_xie

echo Run the various methods on the simulated data.
methods=(sceptre negative_binomial)
theta_sizes=(theta_small theta_correct theta_big)
for i in "${methods[@]}"
do
  for j in "${theta_sizes[@]}"
  do
    Rscript $code_dir"/simulations/run_simulation_generic.R" $offsite_dir $parameter_file_sim $i $j &
  done
done

# Run scMAGeCK separately.
Rscript $code_dir"/simulations/run_simulation_generic.R" $offsite_dir $parameter_file_sim scMAGeCK &
