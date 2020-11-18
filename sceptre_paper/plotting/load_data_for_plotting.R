offsite_dir = "/Users/ekatsevi/Dropbox (Penn)/project-files/sceptre_files"

# Gasperini results
original_results_gasp = read_tsv(sprintf("%s/data/gasperini/raw/GSE120861_all_deg_results.at_scale.txt", offsite_dir))
resampling_results_gasp = read.fst(sprintf("%s/results/gasperini/sceptre/resampling_results.fst", offsite_dir)) %>% as_tibble()
likelihood_results_gasp = read.fst(sprintf("%s/results/gasperini/negative_binomial/likelihood_results.fst", offsite_dir)) %>% as_tibble()
covariates_gasp = read.fst(sprintf("%s/data/gasperini/processed/cell_covariate_model_matrix.fst", offsite_dir)) %>% as_tibble()

# Xie results
original_results_xie = readRDS(sprintf("%s/data/xie/processed/xie_p_values.rds", offsite_dir))
resampling_results_xie = read.fst(sprintf("%s/results/xie/sceptre/all_results_with_names.fst", offsite_dir)) %>% as_tibble()
likelihood_results_xie = read.fst(sprintf("%s/results/xie/negative_binomial/all_results.fst", offsite_dir)) %>% as_tibble()
covariates_xie = read.fst(sprintf("%s/data/xie/processed/covariate_model_matrix.fst", offsite_dir)) %>% as_tibble()
grna_indicator_matrix_xie = read.fst(sprintf("%s/data/xie/processed/grna_indicator_matrix.fst", offsite_dir)) %>% as_tibble()

# simulation results
simulation_results = read.fst(sprintf("%s/results/simulations/all_results.fst", offsite_dir)) %>% as_tibble()

# enrichment results
rejected_pairs_HIC = read_tsv(sprintf("%s/results/gasperini/enrichment/rejected_pairs_HIC.tsv", offsite_dir),
                              col_types = "cciiiciilliid")
TF_enrichments = read_tsv(sprintf("%s/results/gasperini/enrichment/TF_enrichments.tsv", offsite_dir), 
                          col_types = "ccd")
paired_fractions = read_tsv(sprintf("%s/results/gasperini/enrichment/TF_paired_enhancer_fractions.tsv", offsite_dir),
                            col_types = "cidd")

# dispersion coefficients alpha_0 and alpha_1 from equation (6) in DESeq2 paper
disp_coeffs = as.numeric(readRDS(sprintf("%s/data/gasperini/processed/disp_coefficients.rds", offsite_dir)))

# dispersion table (mean expressions and dispersion estimates for each gene)
disp_table = read.fst(sprintf("%s/data/gasperini/processed/disp_table.fst", offsite_dir)) %>% as_tibble()
