# offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
small_example <- FALSE

processed_dir <- paste0(offsite_dir, "/data/xie/processed")
results_dir_negbin <- paste0(offsite_dir, "/results/xie/negative_binomial")
bulk_regions <- readRDS(paste0(processed_dir, "/bulk_region_names.rds"))
gRNA_gene_pairs <- read.fst(paste0(processed_dir, "/gRNA_gene_pairs.fst"))
covariate_matrix <- read.fst(paste0(processed_dir, "/covariate_model_matrix.fst"))
cell_gene_expression_matrix_info <- readRDS(paste0(processed_dir, "/exp_mat_sub_metadata.rds")) 
cell_gene_expression_matrix_info$backingfile <- paste0(processed_dir, "/expression_matrix_sub")
cell_gene_expression_matrix <- cell_gene_expression_matrix_info %>% load_fbm
ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
gRNA_indicator_matrix_fp <- paste0(processed_dir, "/gRNA_indicator_matrix.fst")
regularization_amount <- 3
cell_subset <- seq(1, nrow(cell_gene_expression_matrix))
seed <- 1234
B <- 500
pod_sizes <- c(gene = 100, gRNA = 100, pair = 500)
storage_location <- c(gene_precomp_dir = paste0(offsite_dir, "/data/xie/precomp/gene"), gRNA_precomp_dir = paste0(offsite_dir, "/data/xie/precomp/gRNA"), results_dir = paste0(offsite_dir, "/results/xie/sceptre"), log_dir = paste0(offsite_dir, "/logs/xie"))
gene_precomp_dir <- storage_location[["gene_precomp_dir"]]
gRNA_precomp_dir <- storage_location[["gRNA_precomp_dir"]]
results_dir <- storage_location[["results_dir"]]
log_dir <- storage_location[["log_dir"]]

if (small_example) {
  gRNA_gene_pairs <- slice_sample(gRNA_gene_pairs, n = 20)
  pod_sizes <- c(gene = 2, gRNA = 5, pair = 5)
}


# test ARL15-enh against ARL15
if (FALSE) {
  expressions <- cell_gene_expression_matrix[, ordered_gene_ids == "ENSG00000185305.10"]
  indicators <- read_fst(gRNA_indicator_matrix_fp, columns = "chr5:54325645-54326045") %>% dplyr::pull() %>% as.integer()
  fit <- run_sceptre_gRNA_gene_pair(expressions = expressions,
                             gRNA_indicators = indicators,
                             covariate_matrix = covariate_matrix,
                             B = 500) 
}
