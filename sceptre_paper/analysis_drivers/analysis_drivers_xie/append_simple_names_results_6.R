args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "~/research_code/sceptre-manuscript/" else args[1]
require(fst)
require(tidyverse)
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

p_vals_sceptre <- paste0(results_dir, "/all_results.fst") %>% read.fst()

gene_ids_in_use <- p_vals_sceptre$gene_id %>% as.character()
all_gene_names <- paste0(processed_dir, "/ordered_genes.RDS") %>% readRDS()
all_gene_ids <- paste0(processed_dir, "/ordered_gene_ids.RDS") %>% readRDS()
gene_names <- all_gene_names[match(x = gene_ids_in_use, table = all_gene_ids)]
  
bulk_region_names <- paste0(processed_dir, "/bulk_region_names.rds") %>% readRDS()
enh_names <- bulk_region_names$region_name[match(x = as.character(p_vals_sceptre$gRNA_id), table = bulk_region_names$region)]
p_vals_sceptre_modified <- mutate(p_vals_sceptre, gene_names = factor(gene_names), enh_names = factor(enh_names))

write.fst(x = p_vals_sceptre_modified, path = paste0(results_dir, "/all_results_with_names.fst"))
