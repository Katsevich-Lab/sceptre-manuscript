args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
# Define the code and offsite dirs
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_gasp/sceptre_function_args.R" else args[2]
# source the function arguments
source(param_file)

dictionaries_already_created <- as.logical(args[3])
# Option: are the dictionaries already in place? If so, simply load the dictionaries and print to standard output the number of pods.

if (dictionaries_already_created) {
  n_gene_pods <- paste0(gene_precomp_dir, "/gene_dictionary.fst") %>% read.fst(columns = "pod_id") %>% pull() %>% max()
  n_gRNA_pods <- paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst") %>% read.fst(columns = "pod_id") %>% pull() %>% max()
  n_pair_pods <- paste0(results_dir, "/results_dictionary.fst") %>% read.fst(columns = "pod_id") %>% pull() %>% max()
  paste(n_gene_pods, n_gRNA_pods, n_pair_pods) %>% cat
} else {
  x <- if (!is.null(log_dir)) file.remove(list.files(log_dir, full.names = TRUE))
  dicts <- create_and_store_dictionaries(gRNA_gene_pairs = gRNA_gene_pairs, gene_precomp_dir = gene_precomp_dir, gRNA_precomp_dir = gRNA_precomp_dir, results_dir = results_dir, pod_sizes = pod_sizes)
  paste(dicts$n_pods[["gene"]], dicts$n_pods[["gRNA"]], dicts$n_pods[["pairs"]]) %>% cat
}