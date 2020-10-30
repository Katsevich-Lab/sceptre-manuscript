#' Run negative binomial model (known size)
#'
#' Runs NB model on gene-gRNA pair where gene has known size. Returns the one-sided p-value.
#'
#' @param expressions vector of expressions
#' @param gRNA_indicators vector of gRNA indicators
#' @param covariate_matrix covariate matrix
#' @param gene_size precomputed gene size
#'
#' @return a left-tailed p-value
#' @export
#'
#' @examples
#' offsite_dir <- "/Volumes/tims_new_drive/research/sceptre_files"
#' processed_dir <- paste0(offsite_dir, "/data/xie/processed")
#' covariate_matrix <- read.fst(paste0(processed_dir, "/covariate_model_matrix.fst"))
#' ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
#' cell_gene_expression_matrix <- readRDS(paste0(processed_dir, "/exp_mat_metadata.rds")) %>% load_fbm
#' expressions <- cell_gene_expression_matrix[, which(ordered_gene_ids == "ENSG00000197530.12")]
#' gRNA_indicators <- paste0(processed_dir, "/gRNA_indicator_matrix.fst") %>% read.fst() %>% pull()
#' gene_size <- readRDS(file = paste0(offsite_dir, "/data/xie/precomp/gene/gene_size_unreg_2.rds"))[["ENSG00000197530.12"]]
run_NB_model_known_size <- function(expressions, gRNA_indicators, covariate_matrix, gene_size) {
  full_covariate_matrix <- mutate(covariate_matrix, gRNA_indicator = as.integer(gRNA_indicators))
  fit <- vglm(formula = expressions ~ ., family = negbinomial.size(gene_size), data = full_covariate_matrix)
  z_val <- summaryvglm(fit)@coef3["gRNA_indicator", "z value"]
  p <- pnorm(q = z_val, lower.tail = TRUE)
  return(p)
}


#' Run NB model at scale
#'
#' Runs the negative binomial regression model at scale.
#'
#' @export
run_nb_model_at_scale <- function(pod_id, gene_precomp_dir, results_dir_negbin, cell_gene_expression_matrix, ordered_gene_ids, gRNA_indicator_matrix_fp, covariate_matrix, cell_subset, log_dir) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/result_negbinmom_", pod_id, ".Rout"))
  results_dict <- read.fst(paste0(results_dir_negbin, "/results_dictionary.fst")) %>% filter(pod_id == !!pod_id)
  to_save_fp <- results_dict$result_file[1] %>% as.character()
  gene_ids <- results_dict$gene_id
  gene_sizes <- readRDS(paste0(gene_precomp_dir, "/size_reg_file.rds"))[gene_ids]
  p_vals <- sapply(1:nrow(results_dict), function(i) {
    curr_gene <- results_dict[[i, "gene_id"]] %>% as.character()
    curr_gRNA <- results_dict[[i, "gRNA_id"]] %>% as.character()
    cat(paste0("Running NB regression of gene ", curr_gene, " on gRNA ", curr_gRNA, ".\n"))

    expressions <- cell_gene_expression_matrix[, which(curr_gene == ordered_gene_ids)]
    gRNA_indicators <- read.fst(path = gRNA_indicator_matrix_fp, columns = curr_gRNA) %>% pull() %>% as.integer()
    gene_size <- gene_sizes[[curr_gene]]

    if (!is.null(cell_subset)) {
      expressions <- expressions[cell_subset]
      gRNA_indicators <- gRNA_indicators[cell_subset]
      covariate_matrix <- covariate_matrix[cell_subset,]
    }
    run_NB_model_known_size(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = covariate_matrix, gene_size = gene_size)
  })
  out <- select(results_dict, gRNA_id, gene_id) %>% mutate(p_value = p_vals)
  write_fst(x = out, path = to_save_fp)
}
