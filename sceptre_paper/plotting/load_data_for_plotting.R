offsite_dir = "/Users/ekatsevi/Dropbox (Penn)/project-files/SCEPTRE"

original_results = read_tsv(sprintf("%s/data/raw/CRISPR/GSE120861_all_deg_results.at_scale.txt", offsite_dir))
resampling_results = read.fst(sprintf("%s/resampling_results.fst", offsite_dir)) %>% as_tibble()
likelihood_results = read.fst(sprintf("%s/likelihood_results.fst", offsite_dir)) %>% as_tibble()

covariates_xie = read.fst(sprintf("%s/cell_covariate_model_matrix_xie.fst", offsite_dir)) %>% as_tibble()
covariates_gasp = read.fst(sprintf("%s/cell_covariate_model_matrix_gasp.fst", offsite_dir)) %>% as_tibble()

grna_indicator_matrix_xie = read.fst(sprintf("%s/grna_indicator_matrix_xie.fst", offsite_dir)) %>% as_tibble()
covariates_xie$guide_count = rowSums(grna_indicator_matrix_xie)
