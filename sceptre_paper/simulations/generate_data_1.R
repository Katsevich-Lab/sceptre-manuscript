args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
simulation_data_dir <- paste0(offsite_dir, "/data/simulations")

# We generate datasets corresponding to 4 settings: no zero inflation with theta = 0.2, 1, 5; and 0.25 zero inflation with theta = 1.
dataset_parameters <- tibble(theta = c(0.2, 1, 5, 1), zero_inflation = c(0, 0, 0, 0.25))

n_sim <- 500
n_cells <- 1000
# We generate n_sim random datasets of n_cells each.
gasp_cov_mat <- paste0(offsite_dir, "/data/gasperini/processed/cell_covariate_model_matrix.fst") %>% read.fst() %>% select(lg_total_umis, prep_batch) %>% mutate(prep_batch = factor(prep_batch, levels = c("prep_batch_1", "prep_batch_2"), labels = c("prep_batch_1", "prep_batch_2")))
set.seed(1234)
# First, generate the covariate matrices iid.
covariate_matrices <- replicate(n = n_sim, expr = {
  gasp_cov_mat %>% sample_n(n_cells, replace = TRUE)
}, simplify = FALSE)
# Next, generate the xs and ys
covariate_effects_gRNA <- c(-7, 0.5, -2)
covariate_effects_gene <- c(-2.5, 0.5, -2)
n_enh <- 5
seed <- 1234

datasets <- map(.x = set_names(1:nrow(dataset_parameters), paste0("dataset_", 1:nrow(dataset_parameters))), .f = function(i) {
  curr_theta <- dataset_parameters[[i,"theta"]]
  curr_zero_inf <- dataset_parameters[[i, "zero_inflation"]]
  curr_sim_data <- simulate_crispr_screen_data(covariate_matrices = covariate_matrices, covariate_effects_gRNA = covariate_effects_gRNA, covariate_effects_gene = covariate_effects_gene, zero_inflation = curr_zero_inf, theta = curr_theta, n_enh = n_enh, seed = seed)
  list(theta = curr_theta, zero_inf = curr_zero_inf, ys = curr_sim_data$ys, xs = curr_sim_data$xs)
})

saveRDS(object = covariate_matrices, file = paste0(simulation_data_dir, "/covariate_matrices.rds"))
saveRDS(object = datasets, file = paste0(simulation_data_dir, "/simulated_data.rds"))
