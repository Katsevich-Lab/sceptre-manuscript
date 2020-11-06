args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_xie/sceptre_function_args.R" else args[2]
source(param_file)

simulation_dir <- paste0(offsite_dir, "/data/simulations")

# Load the selected gene and gRNA
chosen_gene_gRNA <- readRDS(paste0(simulation_dir, "/simulation_gRNA_gene.rds"))
my_gene <- chosen_gene_gRNA[["gene"]]
my_gRNA <- chosen_gene_gRNA[["gRNA"]]
cell_subsets <- readRDS(paste0(processed_dir, "/cell_subsets.rds"))

# obtain the data
all_expressions <- cell_gene_expression_matrix[, which(my_gene == ordered_gene_ids)]
all_indicators <- read.fst(gRNA_indicator_matrix_fp, columns = my_gRNA) %>% pull() %>% as.integer()
simulation_covariate_matrix <- select(covariate_matrix, -log_n_gRNA_umis)

# subset the exploratory data
expressions_exploratory <- all_expressions[cell_subsets$exploratory_cells]
gRNA_indicators_exploratory <- all_indicators[cell_subsets$exploratory_cells]
covariate_matrix_exploratory <- simulation_covariate_matrix[cell_subsets$exploratory_cells,]

# Fit the GLMs for the expression and indicator data. Note that we are under the null hypothesis (i.e., no effect of gRNA on gene expression).
fit_nb_exploratory <- glm.nb(formula = expressions_exploratory ~ ., data = covariate_matrix_exploratory)
fit_logistic_exploratory <- glm(gRNA_indicators_exploratory ~ ., family = binomial(), data = covariate_matrix_exploratory)

# Obtain the fitted means of the nb and logistic models on the confirmatory data. For the nb regression, also obtain theta.
covariate_matrix_validation <- simulation_covariate_matrix[cell_subsets$validation_cells,]
fitted_nb_means <- predict.glm(object = fit_nb_exploratory, newdata = covariate_matrix_validation, type = "response") %>% as.numeric()
fitted_theta <- fit_nb_exploratory$theta
fitted_logistic_probs <- predict.glm(object = fit_logistic_exploratory, newdata = covariate_matrix_validation, type = "response") %>% as.numeric()

# Randomly generate n_sim = 500 datasets, where we have n_enh = 20 enhancers
set.seed(1234)
n_sim <- 5
n_enh <- 5

synthetic_ys <- map_dfc(.x = set_names(1:n_sim, paste0("rep_", 1:n_sim)), .f = function(i) {
  cat(paste0("Generating y vector ", i, ".\n"))
  sapply(X = fitted_nb_means, FUN = function(curr_mean) rnegbin(n = 1, mu = curr_mean, theta = fitted_theta))
})

synthetic_xs <- map_dfc(.x = 1:n_sim, .f = function(i) {
  cat(paste0("Generating x matrix ", i, ".\n"))
  curr_names <- paste0("rep_", i, "-enh_", 1:n_enh)
  map_dfc(.x = set_names(1:n_enh, curr_names), .f = function(j) {
    sapply(X = fitted_logistic_probs, FUN = function(curr_mean) rbernoulli(n = 1, p = curr_mean))
  })
})

# define different thetas (theta big, theta small, theta_correct)
thetas <- c(theta_small = fitted_theta/2, theta_correct = fitted_theta, theta_big = 2 * fitted_theta)

# save the fitted models
saveRDS(object = list(fit_nb_exploratory = fit_nb_exploratory, fit_logistic_exploratory = fit_logistic_exploratory), file = paste0(simulation_dir, "/fitted_models.rds"))

# save the "validation" covariate matrix, expression matrix, and gRNA indicator matrix
write.fst(synthetic_ys, paste0(simulation_dir, "/expression_matrix.fst"))
write.fst(synthetic_xs, paste0(simulation_dir, "/gRNA_indicator_matrix.fst"))
write.fst(covariate_matrix_validation, paste0(simulation_dir, "/covariate_matrix.fst"))

# save the thetas
saveRDS(thetas, paste0(simulation_dir, "/thetas.rds"))
