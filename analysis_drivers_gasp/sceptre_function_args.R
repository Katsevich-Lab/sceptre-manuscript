# sceptre function arguments; these arguments should be defined in terms of "offsite_dir"
# offsite_dir <- "/Volumes/tims_new_drive/research/sceptre_files"
small_example <- TRUE

processed_dir <- paste0(offsite_dir, "/data/gasperini/processed")
results_dir_negbin <- paste0(offsite_dir, "/results/gasperini/negative_binomial")
gRNA_gene_pairs <- read.fst(paste0(processed_dir, "/gene_gRNA_pairs_to_study.fst"))
covariate_matrix <- read.fst(paste0(processed_dir, "/cell_covariate_model_matrix.fst"))
cell_gene_expression_matrix <- readRDS(paste0(processed_dir, "/expression_FBM_metadata.rds")) %>% load_fbm
ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
gRNA_indicator_matrix_fp <- paste0(processed_dir, "/gRNA_indicators.fst")
regularization_amount <- 3
cell_subset <- readRDS(paste0(processed_dir, "/cells_to_keep.rds"))
seed <- 1234
B <- 500
pod_sizes <- c(gene = 200, gRNA = 200, pair = 200)
storage_location <- c(gene_precomp_dir = paste0(offsite_dir, "/data/gasperini/precomp/gene"), gRNA_precomp_dir = paste0(offsite_dir, "/data/gasperini/precomp/gRNA"), results_dir = paste0(offsite_dir, "/results/gasperini/sceptre"), log_dir = paste0(offsite_dir, "/logs/gasperini"))
gene_precomp_dir <- storage_location[["gene_precomp_dir"]]
gRNA_precomp_dir <- storage_location[["gRNA_precomp_dir"]]
results_dir <- storage_location[["results_dir"]]
log_dir <- storage_location[["log_dir"]]

if (small_example) {
  pod_sizes <- c(gene = 10, gRNA = 10, pair = 10)
  set.seed(1)
  gRNA_gene_pairs <- slice(gRNA_gene_pairs, sample(x = 1:nrow(gRNA_gene_pairs), size = 60, replace = FALSE))
}
