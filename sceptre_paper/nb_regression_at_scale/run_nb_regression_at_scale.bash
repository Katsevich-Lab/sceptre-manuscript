#!/bin/bash

# This script runs the "at-scale" negative binomial regression analysis. Like the at-scale sceptre bash script, this script takes four arguments: (i) the location of the directory containing the run_nb_regression_at_scale.bash script, (ii) the location of the offsite directory, (iii) the name of the R file containing the parameters of the computation, and (iv) the number of processors to use in the computation.

# Note: This bash script is not quite as general as the sceptre bash script. This script serves more of an "in-house" purpose for performing the analyses reported in Katsevich et al 2020. In particular, this script assumes that the sceptre analysis has already been performed, as this script leverages the sceptre precomputation.

nb_at_scale_bash_dir=$1
offsite_dir=$2
parameter_file=$3
n_processors=$4

# Copy the results dictionary from the sceptre results directory into the negative binomial results directory.
n_pods=$(Rscript $nb_at_scale_bash_dir"/create_dictionary.R" $offsite_dir $parameter_file)
wait

# Run the negative binomial regression across all gRNA-gene pairs.
seq 1 $n_pods | xargs -I{} -n 1 -P $n_processors Rscript $nb_at_scale_bash_dir"/run_regression_at_scale.R" $offsite_dir $parameter_file {} &
wait

# Collect the results
Rscript $nb_at_scale_bash_dir"/collect_results.R" $offsite_dir $parameter_file
