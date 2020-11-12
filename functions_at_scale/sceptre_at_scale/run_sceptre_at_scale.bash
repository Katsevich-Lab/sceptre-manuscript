#!/bin/bash

# This script runs an "at-scale" sceptre analysis. Parallelization is done within Unix (instead of within R). The script takes five arguments: (i) the location of the directory containing the run_sceptre_at_scale.bash script, (ii) the location of the offsite directory, (iii) the name of the R file containing the parameters of the computation, (iv) the number of processors to use in the computation, and (v) whether the precomputation already has been completed.

# The parameter_file must define the following variables in the global environment: gene_precomp_dir, gRNA_precomp_dir, results_dir, log_dir, gRNA_gene_pairs, covariate_matrix, cell_gene_expression_matrix, ordered_gene_ids, gRNA_indicator_matrix_fp, cell_subset, seed, B, pod_sizes, select_sizes. These variables can be defined in terms of offsite_dir, which will be available in the global environment.

sceptre_at_scale_bash_dir=$1
offsite_dir=$2
parameter_file=$3
n_processors=$4
precomputation_complete=$5

# Create the precomputation and result dictionaries
pod_sizes=$(Rscript $sceptre_at_scale_bash_dir"/create_dictionaries.R" $offsite_dir $parameter_file $precomputation_complete)
n_gene_pods="$(echo $pod_sizes | cut -d' ' -f1)"
n_gRNA_pods="$(echo $pod_sizes | cut -d' ' -f2)"
n_pair_pods="$(echo $pod_sizes | cut -d' ' -f3)"

if [ $precomputation_complete = "false" ]
then
echo Run the first round of gene precomputations across all gene pods.
seq 1 $n_gene_pods | xargs -I{} -n 1 -P $n_processors Rscript $sceptre_at_scale_bash_dir/"run_gene_precomputation_round_1.R" $offsite_dir $parameter_file {} &
wait

echo Regularize the gene size estimates.
Rscript $sceptre_at_scale_bash_dir"/regularize_gene_size_estimates.R" $offsite_dir $parameter_file
wait

echo Run the second round of gene precomputations across all gene pods.
seq 1 $n_gene_pods | xargs -I{} -n 1 -P $n_processors Rscript $sceptre_at_scale_bash_dir/"run_gene_precomputation_round_2.R" $offsite_dir $parameter_file {} &
wait

echo Run gRNA precomputation across all gRNA pods.
seq 1 $n_gRNA_pods | xargs -I{} -n 1 -P $n_processors Rscript $sceptre_at_scale_bash_dir"/run_gRNA_precomputation.R" $offsite_dir $parameter_file {} &
wait
fi

echo Run gRNA-gene pair analysis across all pair pods.
seq 1 $n_pair_pods | xargs -I{} -n 1 -P $n_processors Rscript $sceptre_at_scale_bash_dir"/run_pair_analysis_at_scale.R" $offsite_dir $parameter_file {} &
wait

echo Collect and save results.
Rscript $sceptre_at_scale_bash_dir"/aggregate_results.R" $offsite_dir $parameter_file
