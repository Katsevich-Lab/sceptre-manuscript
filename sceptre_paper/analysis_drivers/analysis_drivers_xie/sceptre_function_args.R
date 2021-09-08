# offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
small_example <- FALSE
library(sceptre)

processed_dir <- paste0(offsite_dir, "/data/xie/processed")
analysis_ready_dir <- paste0(offsite_dir, "/data/xie/analysis_ready")
results_dir_negbin <- paste0(offsite_dir, "/results/xie/negative_binomial")
bulk_regions <- readRDS(paste0(processed_dir, "/bulk_region_names.rds"))
gRNA_gene_pairs <- fst::read.fst(paste0(processed_dir, "/gRNA_gene_pairs.fst"))
covariate_matrix <- fst::read.fst(paste0(analysis_ready_dir, "/covariate_model_matrix.fst"))
cell_gene_expression_matrix_info <- readRDS(paste0(analysis_ready_dir, "/exp_mat_sub_metadata.rds")) 
cell_gene_expression_matrix_info$backingfile <- paste0(analysis_ready_dir, "/expression_matrix_sub")
cell_gene_expression_matrix <- sceptre::load_fbm(cell_gene_expression_matrix_info)
ordered_gene_ids <- readRDS(paste0(analysis_ready_dir, "/ordered_gene_ids.RDS"))
gRNA_indicator_matrix_fp <- paste0(analysis_ready_dir, "/gRNA_indicator_matrix.fst")
regularization_amount <- 3
cell_subset <- seq(1, nrow(cell_gene_expression_matrix))
seed <- 1234
B <- 500
pod_sizes <- c(gene = 100, gRNA = 100, pair = 500)
storage_location <- c(gene_precomp_dir = paste0(offsite_dir, "/data/xie/precomp/gene"),
                      gRNA_precomp_dir = paste0(offsite_dir, "/data/xie/precomp/gRNA"),
                      results_dir = paste0(offsite_dir, "/results/xie/sceptre"),
                      log_dir = paste0(offsite_dir, "/logs/xie"))
gene_precomp_dir <- storage_location[["gene_precomp_dir"]]
gRNA_precomp_dir <- storage_location[["gRNA_precomp_dir"]]
results_dir <- storage_location[["results_dir"]]
log_dir <- storage_location[["log_dir"]]

# verify dimensions of gRNA matrix, cell gene expression matrix, and covariate matrix are consistent
dim(cell_gene_expression_matrix)
dim(fst::read_fst(gRNA_indicator_matrix_fp))
dim(covariate_matrix)
length(cell_subset)
length(ordered_gene_ids)

if (small_example) {
  gRNA_gene_pairs <- slice_sample(gRNA_gene_pairs, n = 20)
  pod_sizes <- c(gene = 2, gRNA = 5, pair = 5)
}
