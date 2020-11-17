args <- commandArgs(trailingOnly = TRUE)
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
simulation_result_dir <- paste0(offsite_dir, "/results/simulations")
require(sceptre)

fs <- list.files(simulation_result_dir, full.names = TRUE)
res <- fs %>% map(readRDS) %>% reduce(rbind)
write.fst(res, path = paste0(simulation_result_dir, "/all_results.fst"))
