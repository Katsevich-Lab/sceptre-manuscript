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

# Ensure all simulation folders have been created.
Rscript $code_dir"/sceptre_paper/simulations/"check_directories_0.R $code_dir $offsite_dir

echo Generate the simulation data.
Rscript $code_dir"/sceptre_paper/simulations/"generate_data_1.R $offsite_dir

echo Run the various methods on the various simulated datasets.
# methods=(sceptre negative_binomial scMAGeCK)
methods=(scMAGeCK)
datasets=(1 2 3 4)
for i in "${datasets[@]}"
do
  for j in "${methods[@]}"
  do
    Rscript $code_dir"/sceptre_paper/simulations/"run_simulation_generic.R $offsite_dir $i $j &
  done
done
wait

echo Collect the results
Rscript $code_dir"/sceptre_paper/simulations/"collect_results_2.R $offsite_dir
