args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(katsevich2020))
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]

dataset_id <- as.integer(args[2]) # 1
method <- as.character(args[3]) # sceptre

simulation_data_dir <- paste0(offsite_dir, "/data/simulations")
simulation_log_dir <- paste0(offsite_dir, "/logs/simulations")
simulation_result_dir <- paste0(offsite_dir, "/results/simulations")

covariate_matrices <- paste0(simulation_data_dir, "/covariate_matrices.rds") %>% readRDS()
simulated_data <- paste0(simulation_data_dir, "/simulated_data.rds") %>% readRDS() %>% pluck(dataset_id)
f_name <- paste0("sim_res_", method, "_dataset_", dataset_id)
theta <- simulated_data$theta
n_rep <- ncol(simulated_data$ys)
n_enhancers <- ncol(simulated_data$xs)/n_rep

activate_sink(paste0(simulation_log_dir, "/", f_name, ".Rout"))

p_vals <- map_dbl(.x = 1:n_rep, .f = function(i) {
  cat(paste0("Working on replicate ", i, ".\n"))
  # Load the data
  expressions <- simulated_data$ys[[paste0("rep_", i)]]
  if (method %in% c("sceptre", "negative_binomial")) {
    gRNA_indicators <- simulated_data$xs[["rep_1-enh_1"]] %>% as.integer()
  } else if (method == "scMAGeCK") {
    gRNA_indic_matrix <- simulated_data$xs[,paste0("rep_", i, "-enh_", 1:n_enhancers)] %>% apply(., 2, as.integer)
  }
  curr_covariate_matrix <- covariate_matrices[[i]]

  # Run the method
  if (method == "sceptre") {
    out <- run_sceptre_gRNA_gene_pair(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = curr_covariate_matrix, gene_precomp_size = 1, B = 500, seed = 1234) %>% pull(p_value)
  } else if (method == "negative_binomial") {
    out <- run_NB_model_known_size(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = curr_covariate_matrix, gene_size = 1)
  } else if (method == "scMAGeCK") {
    lib_sizes <- as.integer(exp(curr_covariate_matrix$lg_total_umis))
    out <- run_scmageck_simulation(expressions = expressions, gRNA_indic_matrix = gRNA_indic_matrix, lib_sizes = lib_sizes)
  }
  return(out)
})

to_save <- tibble(p_value = p_vals, method = factor(method), theta = factor(theta), zero_inf = factor(simulated_data$zero_inf), dataset_id = factor(dataset_id))
saveRDS(to_save, paste0(simulation_result_dir, "/", f_name, ".rds"))

deactivate_sink()
