args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_gasp/sceptre_function_args.R" else args[2]
source(param_file)

# Determine the theta across genes
files_to_load <- paste0(gene_precomp_dir, "/", grep(pattern = "gene_size_unreg*", x = list.files(gene_precomp_dir), value = TRUE))
all_thetas_unreg <- map(files_to_load, readRDS) %>% reduce(c)
all_thetas_reg <- readRDS(file = paste0(gene_precomp_dir, "/size_reg_file.rds"))

median_gene <- names(which.min(abs(all_thetas_unreg - median(all_thetas))))

expressions <- cell_gene_expression_matrix[, which(ordered_gene_ids == median_gene)]
