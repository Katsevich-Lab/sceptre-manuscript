# Depricated

#' Run gene precomputation at scale (depricated)
#'
#' @param pod_id ID of the gene pod on which to run the precomputation
#' @param gene_precomp_dir file path to the gene precomputation directory
#' @param cell_gene_expression_matrix an FBM containing the expression data (rows cells, columns genes)
#' @param ordered_gene_ids the gene IDs (i.e., names)
#' @param covariate_matrix the cell-specific covariate matrix
#' @param cell_subset (optional) integer vector identifying the cells to use in the model
#' @param gene_sizes (optional) a vector of already-estimated (or known) gene sizes
#' @param log_dir file path to the directory in which to sink the log output
#'
#' @return NULL
#'
#' @examples
#' pod_id <- 1
#' gene_precomp_dir <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/precomp/gene/"
#' cell_gene_expression_matrix_metadata <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/expression_FBM_metadata.rds")
#' cell_gene_expression_matrix <- load_fbm(cell_gene_expression_matrix_metadata)
#' covariate_matrix <- read.fst("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cell_covariate_model_matrix.fst")
#' cell_subset <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cells_to_keep.rds")
#' ordered_gene_ids <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/ordered_gene_ids.RDS")
run_gene_precomputation_at_scale <- function(pod_id, gene_precomp_dir, cell_gene_expression_matrix, ordered_gene_ids, covariate_matrix, cell_subset = NULL, log_dir = NULL, gene_sizes = NULL) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/gene_precomp_", pod_id, ".Rout"))
  # subset covariate matrix by rows
  if (!is.null(cell_subset)) covariate_matrix <- covariate_matrix[cell_subset,]
  # obtain the genes on which to preform the precomputation
  gene_dictionary <- read.fst(paste0(gene_precomp_dir, "/gene_dictionary.fst")) %>% filter(pod_id == !!pod_id)
  gene_ids <- gene_dictionary %>% pull(id) %>% as.character()

  # Run the precomputations
  precomps <- map(gene_ids, function(gene_id) {
    cat(paste0("Running precomputation for gene ", gene_id, ".\n"))
    integer_id <- which(gene_id == ordered_gene_ids)
    expressions <- cell_gene_expression_matrix[,integer_id]
    if (!is.null(cell_subset)) expressions <- expressions[cell_subset]
    curr_gene_size <- gene_sizes[[gene_id]]
    run_gene_precomputation(expressions = expressions, covariate_matrix = covariate_matrix, gene_precomp_size = curr_gene_size)
  })
  names(precomps) <- gene_ids
  out_offsets <- map_dfc(precomps, function(l) l$gene_precomp_offsets)
  out_sizes <- map_dbl(precomps, function(l) l$gene_precomp_size)

  # save the precomputations
  offset_save_fp <- (gene_dictionary %>% pull(offset_file))[1] %>% as.character()
  size_save_fp <- (gene_dictionary %>% pull(size_file))[1] %>% as.character()
  write.fst(out_offsets, offset_save_fp)
  saveRDS(out_sizes, size_save_fp)
  if (!is.null(log_dir)) deactivate_sink()
}

#' Run sceptre using precomputations for gRNAs and genes (old version)
#'
#' This function is the workhorse function of the sceptre package. It runs a distilled CRT using a negative binomial test statistic based on an expression vector, a gRNA indicator vector, an offset vector (from the distillation step), gRNA conditional probabilities, an estimate of the negative binomial size parameter, and the number of resampling replicates.
#'
#' This currently is a one-tailed, left-sided test. Thus, it is suitable for up-regulatory elements like enhancers and promoters but not down-regulatory elements like silencers.
#'
#' @param expressions a vector of gene expressions (in UMI counts)
#' @param gRNA_indicators a vector of gRNA indicators
#' @param gRNA_precomp a vector of conditional probabilities for gRNA assignments
#' @param gene_precomp_size the pre-computed size
#' @param gene_precomp_offsets the pre-computed distillation offsets
#' @param B the number of resamples to make (default 500)
#' @param seed an arguement to set.seed; if null, no seed is set
#'
#' @return a p-value of the null hypothesis of no gRNA effect on gene expression
run_sceptre_using_precomp_old <- function(expressions, gRNA_indicators, gRNA_precomp, gene_precomp_size, gene_precomp_offsets, B, seed, reduced_output = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  # compute the test statistic on the real data
  fit_star <- vglm(formula = expressions[gRNA_indicators == 1] ~ 1, family = negbinomial.size(gene_precomp_size), offset = gene_precomp_offsets[gRNA_indicators == 1])
  t_star <- summaryvglm(fit_star)@coef3["(Intercept)", "z value"]

  # resample B times
  t_nulls <- sapply(1:B, function(i) {
    if (i %% 100 == 0) cat(paste0("Running resample ", i ,"/", B, ".\n"))
    gRNA_indicators_null <- rbinom(n = length(gRNA_precomp), size = 1, prob = gRNA_precomp)
    tryCatch({
      fit_null <- vglm(formula = expressions[gRNA_indicators_null == 1] ~ 1, family = negbinomial.size(gene_precomp_size), offset = gene_precomp_offsets[gRNA_indicators_null == 1])
      summaryvglm(fit_null)@coef3["(Intercept)", "z value"]},
      error = function(e) return(NA),
      warning = function(w) return(NA)
    )
  })
  t_nulls <- t_nulls[!is.na(t_nulls)]

  # Fit a skew-t distribution and obtain a p-value
  p_value_skew_t <- NA
  skew_t_fit <- tryCatch(selm(t_nulls ~ 1, family = "ST"), error = function(e) return(NA), warning = function(w) return(NA))
  if (class(skew_t_fit) == "selm") { # If the fit worked,
    dp <- skew_t_fit@param$dp # then extract the parameters.
    if (!any(is.na(dp))) { # If all the fitted parameters are numbers,
      p_value_skew_t <- pst(x = t_star, dp = dp) # then compute the skew t-based p-value.
    }
  }
  # check if the skew-t fit worked
  skew_t_fit_success <- !is.na(p_value_skew_t)
  if (skew_t_fit_success) {
    out_p <- p_value_skew_t
    skew_t_mle <- dp
  } else {
    out_p <- mean(c(-Inf, t_nulls) <= t_star)
    skew_t_mle <- c(xi = NA, omega = NA, alpha = NA, nu = NA)
  }

  # Prepare the output
  if (reduced_output) {
    out <- data.frame(p_value =  out_p, skew_t_fit_success = skew_t_fit_success, xi = skew_t_mle[["xi"]], omega = skew_t_mle[["omega"]], alpha = skew_t_mle[["alpha"]], nu = skew_t_mle[["nu"]])
  } else {
    out <- list(p_value = out_p,
                skew_t_fit_success = skew_t_fit_success,
                skew_t_mle = skew_t_mle,
                z_value = t_star,
                resampled_z_values = t_nulls)
  }
  return(out)
}

