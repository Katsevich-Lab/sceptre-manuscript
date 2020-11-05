args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(katsevich2020))
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
simulation_param_file <- if (is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/simulations/simulation_parameters.R" else args[2]
source(simulation_param_file)

method <- args[3] # "sceptre"
theta_size <- args[4] # "theta_small"

theta <- if (!is.na(theta_size)) thetas[[theta_size]] else NA
n_rep <- ncol(read.fst(path = expression_matrix_fp, from = 1, to = 1))
n_enhancers <- ncol(read.fst(path = gRNA_indicator_matrix_fp, from = 1, to = 1))/n_rep

p_vals <- map_dbl(.x = 1:n_rep, .f = function(i) {
  cat(paste0("Working on replicate ", i, ".\n"))
  # Load the data
  expressions <- read.fst(expression_matrix_fp, paste0("rep_", i)) %>% pull()
  if (method %in% c("sceptre", "negative_binomial")) {
    gRNA_indicators <- read.fst(gRNA_indicator_matrix_fp, paste0("rep_", i, "-enh_1")) %>% pull() %>% as.integer()
  } else if (method == "scMAGeCK") {
    gRNA_indic_matrix <- read.fst(gRNA_indicator_matrix_fp, paste0("rep_", i, "-enh_", 1:n_enhancers)) %>% apply(., 2, as.integer)
  }

  # Run the method
  if (method == "sceptre") {
    out <- run_sceptre_gRNA_gene_pair(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = covariate_matrix, gene_precomp_size = theta, B = B, seed = seed) %>% pull(p_val)
  } else if (method == "negative_binomial") {
    out <- run_NB_model_known_size(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = covariate_matrix, gene_size = theta)
  } else if (method == "scMAGeCK") {
    lib_sizes <- as.integer(exp(covariate_matrix$log_n_umis))
    out <- run_scmageck_simulation(expressions = expressions, gRNA_indic_matrix = gRNA_indic_matrix, lib_sizes = lib_sizes)
  }
})


to_save <- tibble(p_value = p_vals, method = factor(method), theta = factor(theta), theta_size = factor(theta_size))
f_name <- paste0("sim_res_", if (!is.na(theta_size)) paste0(theta_size,"_") else "", method, ".rds")
saveRDS(to_save, paste0(simulation_dir, "/", f_name))
