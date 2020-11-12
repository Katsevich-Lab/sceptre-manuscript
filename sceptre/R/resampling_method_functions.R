#' Fit skew-t
#'
#' Fits a skew-t distribution on a set of resampled test statistics.
#' @param t_nulls
#'
#' @return
#' A list containing (i) skew_t_fit_success (boolean), (ii) out_p (the p-value), and (iii) skew-t mle (a vector containing the fitted MLEs, NA if fit failed).
#' @export

fit_skew_t <- function(t_nulls, t_star) {
  p_value_skew_t <- NA
  skew_t_fit <- tryCatch(selm(t_nulls ~ 1, family = "ST"), error = function(e) return(NA), warning = function(w) return(NA))
  if (class(skew_t_fit) == "selm") { # If the fit worked,
    dp <- skew_t_fit@param$dp # then extract the parameters.
    if (!any(is.na(dp))) { # If all the fitted parameters are numbers,
      p_value_skew_t <- pmax(.Machine$double.eps, pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)) # then compute the skew t-based p-value. pst(x = t_star, dp = dp)
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
  return(list(skew_t_fit_success = skew_t_fit_success, out_p = out_p, skew_t_mle = skew_t_mle))
}

#' Run sceptre using precomputations for gRNAs and genes.
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
#' @export
#' @return a p-value of the null hypothesis of no gRNA effect on gene expression
#' @examples
#' offsite_dir <- "/Volumes/tims_new_drive/research/sceptre_files"
#' source("/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_gasp/sceptre_function_args.R")
#' gene_id <- "ENSG00000008256"
#' gRNA_id <- "ACTB_TSS"
#' expressions <- cell_gene_expression_matrix[, which(ordered_gene_ids == gene_id)][cell_subset]
#' gRNA_indicators <- (read.fst(gRNA_indicator_matrix_fp, columns = "ASAH1_TSS") %>% pull() %>% as.integer())[cell_subset]
#' gene_sizes <- readRDS(paste0(gene_precomp_dir, "/size_reg_file.rds"))
#' gene_precomp_size <- gene_sizes[[gene_id]]
#' gRNA_precomp <- paste0(gRNA_precomp_dir, "/gRNA_precomp_1.fst") %>% read.fst(columns = gRNA_id) %>% pull()
#' gene_precomp_offsets <- paste0(gene_precomp_dir, "/gene_offsets_1.fst") %>% read.fst(columns = gene_id) %>% pull()
#' run_sceptre_using_precomp(expressions, gRNA_indicators, gRNA_precomp, gene_precomp_size, gene_precomp_offsets, 500, 1234)
run_sceptre_using_precomp <- function(expressions, gRNA_indicators, gRNA_precomp, gene_precomp_size, gene_precomp_offsets, B, seed, reduced_output = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  # compute the test statistic on the real data
  fit_star <- vglm(formula = expressions[gRNA_indicators == 1] ~ 1, family = negbinomial.size(gene_precomp_size), offset = gene_precomp_offsets[gRNA_indicators == 1])
  t_star <- summaryvglm(fit_star)@coef3["(Intercept)", "z value"]

  # Define a closure to resample B times (omitting the NAs)
  resample_B_times <- function(my_B) {
    t_nulls <- sapply(1:my_B, function(i) {
      if (i %% 100 == 0) cat(paste0("Running resample ", i ,"/", my_B, ".\n"))
      gRNA_indicators_null <- rbinom(n = length(gRNA_precomp), size = 1, prob = gRNA_precomp)
      tryCatch({
        fit_null <- vglm(formula = expressions[gRNA_indicators_null == 1] ~ 1, family = negbinomial.size(gene_precomp_size), offset = gene_precomp_offsets[gRNA_indicators_null == 1])
        summaryvglm(fit_null)@coef3["(Intercept)", "z value"]},
        error = function(e) return(NA),
        warning = function(w) return(NA)
      )
    })
    t_nulls[!is.na(t_nulls)]
  }

  # resample B times
  t_nulls <- resample_B_times(B)

  # obtain a p-value
  skew_t_fit <- fit_skew_t(t_nulls, t_star)

  # determine if the fit was successful
  if (!skew_t_fit$skew_t_fit_success || length(t_nulls) <= floor(0.95 * B)) { # If the skew-t fit failed, then try again.
    t_nulls_second_set <- resample_B_times(9 * B)
    t_nulls <- c(t_nulls, t_nulls_second_set)
    skew_t_fit <- fit_skew_t(t_nulls, t_star)
  }

  # Prepare the output
  if (reduced_output) {
    out <- data.frame(p_value = skew_t_fit$out_p, skew_t_fit_success = skew_t_fit$skew_t_fit_success, xi = skew_t_fit$skew_t_mle[["xi"]], omega = skew_t_fit$skew_t_mle[["omega"]], alpha = skew_t_fit$skew_t_mle[["alpha"]], nu = skew_t_fit$skew_t_mle[["nu"]], z_value = t_star, n_successful_resamples = length(t_nulls))
  } else {
    out <- list(p_value = skew_t_fit$out_p,
                skew_t_fit_success = skew_t_fit$skew_t_fit_success,
                skew_t_mle = skew_t_fit$skew_t_mle,
                z_value = t_star,
                resampled_z_values = t_nulls)
  }
  return(out)
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
#' @export
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


#' Run sceptre on a gRNA-gene pair
#'
#' This function runs the sceptre algorithm on a single gRNA-gene pair. It requires as arguments the gene expression vector, the gRNA indicator vector, and the covariate matrix. Users optionally can pass the gRNA precomputation or gene precomputation as arguments.
#'
#' @param expressions a vector a gene expressions
#' @param gRNA_indicators a vector of gRNA inicators
#' @param covariate_matrix the matrix of cell-specific covariates (e.g., library size, batch effect, cell cycle, etc.)
#' @param gene_precomp_size (optional) the pre-computed size of the gene NB distribution
#' @param B number of resamples (default 500)
#' @param seed (optional) seed to the random number generator
#'
#' @return a p-value of the null hypothesis of no gRNA effect on gene expression
#' @export
#'
#' @examples
#' offsite_dir <- "/Volumes/tims_new_drive/research/sceptre_files"
#' source("/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_xie/sceptre_function_args.R")
#' expressions <- cell_gene_expression_matrix[,2][cell_subset]
#' gRNA_indicators <- (read.fst(gRNA_indicator_matrix_fp, columns = "chr5:54325645-54326045") %>% pull() %>% as.integer())[cell_subset]
#' covariate_matrix <- if (nrow(covariate_matrix) == 106666) covariate_matrix else covariate_matrix[cell_subset,]
#' run_sceptre_gRNA_gene_pair(expressions, gRNA_indicators, covariate_matrix)
#'
run_sceptre_gRNA_gene_pair <- function(expressions, gRNA_indicators, covariate_matrix, gene_precomp_size = NULL, B = 500, seed = NULL) {
  cat(paste0("Running gRNA precomputation.\n"))
  gRNA_precomp <- run_gRNA_precomputation(gRNA_indicators, covariate_matrix)

  cat(paste0("Running gene precomputation.\n"))
  gene_precomp <- run_gene_precomputation(expressions, covariate_matrix, gene_precomp_size)

  out <- run_sceptre_using_precomp(expressions = expressions, gRNA_indicators = gRNA_indicators, gRNA_precomp = gRNA_precomp, gene_precomp_size = gene_precomp$gene_precomp_size, gene_precomp_offsets = gene_precomp$gene_precomp_offsets, B = B, seed = seed)
  return(out)
}


#' Run gRNA precomputation
#'
#' This function runs the precomputation for a given gRNA.
#'
#' @param gRNA_indicators a vector of gRNA indicators
#' @param covariate_matrix the cell-specific covariate matrix
#'
#' @export
#' @return the fitted probabilities
run_gRNA_precomputation <- function(gRNA_indicators, covariate_matrix) {
  fit_model_grna <- glm(gRNA_indicators ~ ., family = binomial(), data = covariate_matrix)
  out <- as.numeric(fitted(fit_model_grna))
  return(out)
}


#' Run gene precomputation
#'
#' This function runs the precomputation for a given gene. In particlar, it fits an NB regression of expression against covariates. The estimate of theta (i.e., the NB size parameter) is obtained from glm.nb function. This is sensible as, under the null hypothesis, the NB model without the gRNA indicator is true. Offsets are obtained by log-transforming the fitted values.
#'
#' @param expressions the vector of gene expressions
#' @param covariate_matrix the cell-specific covariate matrix
#' @param gene_precomp_size the pre-computed size parameter (NULL if none)
#'
#' @return a named list containing two items: offsets and size.
#' @export
run_gene_precomputation <- function(expressions, covariate_matrix, gene_precomp_size) {
  # cases on gene_precomp_size
  if (is.null(gene_precomp_size)) { # no size supplied; use glm.nb to estimate size and fit model

    backup_2 <- function(pois_fit) {
      theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
    }
    backup <- function() {
      pois_fit <- glm(expressions ~ ., data = covariate_matrix, family = poisson())
      gene_precomp_size_out <- tryCatch({
        theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1]
      }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
      fit_nb <- vglm(formula = expressions ~ ., family = negbinomial.size(gene_precomp_size_out), data = covariate_matrix)
      fitted_vals <- as.numeric(fittedvlm(fit_nb))
      list(fitted_vals = fitted_vals, gene_precomp_size_out = gene_precomp_size_out)
    }

    result <- tryCatch({
      fit_nb <- glm.nb(formula = expressions ~ ., data = covariate_matrix)
      fitted_vals <- as.numeric(fit_nb$fitted.values)
      gene_precomp_size_out <- fit_nb$theta
      list(fitted_vals = fitted_vals, gene_precomp_size_out = gene_precomp_size_out)
    }, error = function(e) backup(), warning = function(w) backup())

    fitted_vals <- result$fitted_vals; gene_precomp_size_out <- result$gene_precomp_size_out

  } else { # size supplied; use vglm to fit model
    gene_precomp_size_out <- gene_precomp_size
    fit_nb <- vglm(formula = expressions ~ ., family = negbinomial.size(gene_precomp_size_out), data = covariate_matrix)
    fitted_vals <- as.numeric(fittedvlm(fit_nb))
  }

  gene_precomp_offsets_out <- log(fitted_vals)
  out <- list(gene_precomp_offsets = gene_precomp_offsets_out, gene_precomp_size = gene_precomp_size_out)
  return(out)
}
