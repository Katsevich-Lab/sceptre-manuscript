#' Simulate CRISPR screen data
#'
#' @param covariate_matrices a list of covariate matrices
#' @param covariate_effects_gRNA the effect of the technical covariates on gRNA presence probability (a vector of values, whose length corresponds to the number of covariates; must include an intercept term)
#' @param covariate_effects_gene the effect of the technical covariates on mean gene expression (a vector of values, whose length corresponds to the number of covariates; must include an intercept term)
#' @param zero_inflation zero inflation rate (between 0 and 1)
#' @param theta size of the negative binomial distribution for gene expression
#' @param n_enhancers number of enhancers to sample (conditionally on the covariates)
#' @param seed (optonal) seed to the random number generator
#'
#' @return list with components (i) cell_by_gene_expressions, (ii) cell_by_enhancer_perturbation_indicators, and (iii) Z
#' @export
#'
#' @examples
#' n_cells <- 1000
#' gasp_cov_mat <- read.fst("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cell_covariate_model_matrix.fst")
#' covariate_matrices <- replicate(n = 500, expr = {
#'  gasp_cov_mat %>% select(lg_total_umis, prep_batch) %>% mutate(prep_batch = factor(prep_batch, levels = c("prep_batch_1", "prep_batch_2"), labels = c("prep_batch_1", "prep_batch_2"))) %>% sample_n(n_cells)
#' }, simplify = FALSE)
#' mean_log_umis <- mean(gasp_cov_mat$lg_total_umis)
#' mean_batch <- mean(ifelse(gasp_cov_mat$prep_batch == "prep_batch_1", 0, 1))
#' covariate_effects_gRNA <- c(-7, 0.5, -2)
#' covariate_effects_gene <- c(-2.5, 0.5, -2)
#' mean_mu_gRNA <- (covariate_effects_gRNA[1] + covariate_effects_gRNA[2] * mean_log_umis + covariate_effects_gRNA[3] * mean_batch) %>% binomial()$linkinv()
#' mean_mu_gene <- (covariate_effects_gene[1] + covariate_effects_gene[2] * mean_log_umis + covariate_effects_gene[3] * mean_batch) %>% negative.binomial(1)$linkinv()
#' zero_inflation <- 0
#' theta <- 1
#' n_enh <- 5
#' seed <- 1234
#' sim_data <- simulate_crispr_screen_data(covariate_matrices, covariate_effects_gRNA, covariate_effects_gene, zero_inflation, theta, n_enh, seed)
simulate_crispr_screen_data <- function(covariate_matrices, covariate_effects_gRNA, covariate_effects_gene, zero_inflation, theta, n_enh, seed) {
  if (!is.null(seed)) set.seed(seed)
  # First, generate a model matrix for each input covariate matrix
  model_matrices <- map(.x = covariate_matrices, .f = function(covariate_matrix) {
    vars <- colnames(covariate_matrix)
    form <- paste0("~ ", paste0(vars, collapse = " + ")) %>% as.formula()
    model_mat <- model.matrix(form, covariate_matrix)
  })
  # define a quick helper function to compute GLM means as a function of coefficients and model matrices
  get_glm_means <- function(model_matrices, coefs, link_inv) {
    map(.x = model_matrices, function(model_matrix) {
      ls <- (model_matrix %*% coefs) %>% as.numeric()
      link_inv(ls)
    })
  }
  # Get the means of the genes and gRNAs
  gRNA_means <- get_glm_means(model_matrices, covariate_effects_gRNA, link_inv = binomial()$linkinv)
  gene_means <- get_glm_means(model_matrices, covariate_effects_gene, link_inv = negative.binomial(theta = theta)$linkinv)
  n_sim <- length(covariate_matrices)

  # Create random samples for the genes based on the calculated means
  synthetic_ys <- map_dfc(.x = set_names(1:n_sim, paste0("rep_", 1:n_sim)), .f = function(i) {
    cat(paste0("Generating y vector ", i, ".\n"))
    sapply(X = gene_means[[i]], FUN = function(curr_mean) rnegbin(n = 1, mu = curr_mean, theta = theta))
  })

  # Apply the zero-inflation
  indics <- rbernoulli(n = nrow(synthetic_ys) * ncol(synthetic_ys), p = 1 - zero_inflation) %>% as.integer()
  synthetic_ys <- synthetic_ys * matrix(data = indics, nrow = nrow(synthetic_ys), ncol = ncol(synthetic_ys))

  # Likewise, create random samples for the gRNAs based on the calculated means
  synthetic_xs <- map_dfc(.x = 1:n_sim, .f = function(i) {
    cat(paste0("Generating x matrix ", i, ".\n"))
    curr_names <- paste0("rep_", i, "-enh_", 1:n_enh)
    map_dfc(.x = set_names(1:n_enh, curr_names), .f = function(j) {
      sapply(X = gRNA_means[[i]], FUN = function(curr_mean) rbernoulli(n = 1, p = curr_mean))
    })
  })
  return(list(ys = synthetic_ys, xs = synthetic_xs))
}
