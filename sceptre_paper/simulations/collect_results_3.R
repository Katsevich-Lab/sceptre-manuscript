args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
suppressPackageStartupMessages(library(katsevich2020))
suppressPackageStartupMessages(library(scales))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
simulation_param_file <- if (is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/simulations/simulation_parameters.R" else args[2]
source(simulation_param_file)

fs <- list.files(simulation_result_dir, full.names = TRUE)
res <- fs %>% map(readRDS) %>% reduce(rbind)
res <- res %>% mutate(facet_title =  ifelse(is.na(theta_size), "scMAGeCK", as.character(theta_size)) %>% factor(levels = c("theta_small", "theta_correct", "theta_big", "scMAGeCK"), labels = c("Theta too small", "Correct model", "Theta too big", "scMAGeCK")), method_plot = factor(as.character(method), levels = c("negative_binomial", "sceptre", "scMAGeCK"), labels = c("Negative binomial", "SCEPTRE", "scMAGeCK")))
write.fst(res, path = paste0(simulation_result_dir, "/all_results.fst"))
