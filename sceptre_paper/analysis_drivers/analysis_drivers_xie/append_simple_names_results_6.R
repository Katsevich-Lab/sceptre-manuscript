code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

# load data
all_gene_names <- paste0(processed_dir, "/ordered_genes.RDS") %>% readRDS()
all_gene_ids <- paste0(processed_dir, "/ordered_gene_ids.RDS") %>% readRDS()
bulk_region_names <- paste0(processed_dir, "/bulk_region_names.rds") %>% readRDS()
gRNA_gene_pair <- read.fst(paste0(processed_dir, '/gRNA_gene_pairs.fst')) %>% mutate(pair_str = paste0(gRNA_id, "+", gene_id)) # 0 duplicated rows

# function
annotate_results <- function(results) {
  gene_ids_in_use <- results$gene_id %>% as.character()
  gene_names <- all_gene_names[match(x = gene_ids_in_use, table = all_gene_ids)]
  enh_names <- bulk_region_names$region_name[match(x = as.character(results$gRNA_id), table = bulk_region_names$region)]
  all_results <- mutate(results, gene_names = factor(gene_names), enh_names = factor(enh_names))
  resampling_results_xie <- all_results %>% mutate(pair_str = paste0(gRNA_id, "+", gene_id)) # 508 duplicate rows
  idxs_to_remove <- resampling_results_xie %>% select(gene_id, gRNA_id) %>% duplicated() %>% which()
  resampling_results_xie <- resampling_results_xie[-idxs_to_remove,]
  # next, use join to find the result to each gRNA-gene pair
  results_with_type <- left_join(x = select(gRNA_gene_pair, type, pair_str), y = resampling_results_xie, by = "pair_str")
  # finally, apply the fdr correction
  all_results_annotated <- results_with_type %>% group_by(type) %>%
    mutate(adjusted_pvalue = ifelse(type == 'cis', p.adjust(p_value, 'fdr'), NA),
           rejected = ifelse(is.na(adjusted_pvalue), NA, adjusted_pvalue <= 0.1)) %>% ungroup()
  return(all_results_annotated)
}

# annotate three sets of results: sceptre, inb, and monocle nb
p_vals_sceptre <- paste0(results_dir, "/all_results.fst") %>% read.fst()
p_vals_inb <- paste0(results_dir_negative_binomial, "/all_results.fst") %>% read.fst()
p_vals_monocle <- paste0(results_dir_negative_binomial, "/monocle_nb_results.rds") %>% readRDS() %>% rename("p_value" = "p_val")

# annotate the results
sceptre_results_annotated <- annotate_results(p_vals_sceptre)
inb_results_annotated <- annotate_results(p_vals_inb)
monocle_results_annotated <- annotate_results(p_vals_monocle)

# save the annotated results
write_fst(sceptre_results_annotated, paste0(results_dir, "/all_results_annotated.fst"))
write_fst(inb_results_annotated, paste0(results_dir_negative_binomial, "/all_results_annotated.fst"))
write_fst(monocle_results_annotated, paste0(results_dir_negative_binomial, "/monocle_results_annotated.fst"))
