code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
analysis_ready_dir <- paste0(offsite_dir, "data/xie/analysis_ready")

# source args
# Load the expression FBM to compute the gene expressions
exp_mat_metadata <- readRDS(paste0(analysis_ready_dir, "/exp_mat_sub_metadata.rds"))
ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
# update backing file
exp_mat_metadata$backingfile <- paste0(analysis_ready_dir, "/expression_matrix_sub")
# load FBM given backing file
load_fbm <- function(fbm_metadata) {
  out <- bigstatsr::FBM(nrow = fbm_metadata$nrow, ncol = fbm_metadata$ncol,
                        type = fbm_metadata$type, backingfile = fbm_metadata$backingfile,
                        is_read_only = TRUE, create_bk = FALSE)
  return(out)
}
exp_mat <- load_fbm(exp_mat_metadata)
n_cells <- nrow(exp_mat)

# compute the count of each gene
gene_expressions <- bigstatsr::big_apply(X = exp_mat, a.FUN = function(x, ind) {
  colSums(x[,ind])
}, a.combine = "c")

# Finally, put the gene expression counts into a data frame
gene_expressions <- data.frame(gene_id = ordered_gene_ids, expression = gene_expressions)
# saveRDS(gene_expressions, paste0(processed_dir, "/gene_expressions.rds"))

# conduct the binning analysis with 3 bins (fewer pairs were tested)
gene_expressions <- readRDS(paste0(processed_dir, "/gene_expressions.rds"))
results <- fst::read_fst(paste0(results_dir, "/all_results_annotated.fst"))
results_inb <- fst::read_fst(paste0(results_dir_negative_binomial, "/all_results_annotated.fst"))

conduct_binning_analysis <- function(gene_expressions, results, k) {
  results_cis <- results %>% dplyr::filter(type == "cis")
  # divide the cis genes into k groups based on expression level
  cis_genes <- results_cis %>% dplyr::pull(gene_id) %>% as.character() %>% unique()
  gene_expressions_cis <- gene_expressions %>% dplyr::filter(gene_id %in% cis_genes) %>%
    dplyr::mutate(rank = dplyr::ntile(expression, n = k)) %>%
    dplyr::mutate(m_expression = expression/n_cells)
  
  # for each rank, compute the fraction of pairs rejected
  gene_expressions_cis %>% dplyr::group_by(rank) %>% dplyr::group_modify(.f = function(tbl, key) {
    results_cis %>% dplyr::filter(gene_id %in% tbl$gene_id) %>%
      dplyr::summarize(p_rejected = 100 * mean(rejected), n_genes = nrow(tbl),
                       m_exp = mean(tbl$m_expression), n_pairs = length(pair_str))
  })  
}

binning_result_xie <- conduct_binning_analysis(gene_expressions, results, 3)
binning_result_xie_inb <- conduct_binning_analysis(gene_expressions, results_inb, 3)
