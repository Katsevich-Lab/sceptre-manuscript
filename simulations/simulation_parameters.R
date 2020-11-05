# offsite_dir <- "/Volumes/tims_new_drive/research/sceptre_files"
simulation_dir <- paste0(offsite_dir, "/data/simulations")
simulation_log_dir <- paste0(offsite_dir, "/logs/simulations")
simulation_result_dir <- paste0(offsite_dir, "/results/simulations")
thetas <- readRDS(paste0(simulation_dir, "/thetas.rds"))
expression_matrix_fp <- paste0(simulation_dir, "/expression_matrix.fst")
gRNA_indicator_matrix_fp <- paste0(simulation_dir, "/gRNA_indicator_matrix.fst")
covariate_matrix <- paste0(simulation_dir, "/covariate_matrix.fst") %>% read.fst()
seed <- 1234
B <- 500
