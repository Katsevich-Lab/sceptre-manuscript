args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "~/research_code/sceptre-manuscript/" else args[1]
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

# first, add additional columns for gene names and enhancer names
p_vals_sceptre <- paste0(results_dir, "/all_results.fst") %>% read.fst()

gene_ids_in_use <- p_vals_sceptre$gene_id %>% as.character()
all_gene_names <- paste0(processed_dir, "/ordered_genes.RDS") %>% readRDS()
all_gene_ids <- paste0(processed_dir, "/ordered_gene_ids.RDS") %>% readRDS()
gene_names <- all_gene_names[match(x = gene_ids_in_use, table = all_gene_ids)]
  
bulk_region_names <- paste0(processed_dir, "/bulk_region_names.rds") %>% readRDS()
enh_names <- bulk_region_names$region_name[match(x = as.character(p_vals_sceptre$gRNA_id), table = bulk_region_names$region)]
all_results <- mutate(p_vals_sceptre, gene_names = factor(gene_names), enh_names = factor(enh_names))

# create cis file for resampling results for Xie data
gRNA_gene_pair <- read.fst(paste0(processed_dir, '/gRNA_gene_pairs.fst')) %>% mutate(pair_str = paste0(gRNA_id, "+", gene_id)) # 0 duplicated rows
resampling_results_xie <- all_results %>% mutate(pair_str = paste0(gRNA_id, "+", gene_id)) # 508 duplicate rows
# first, remove the duplicated rows from resampling_results_xie
idxs_to_remove <- resampling_results_xie %>% select(gene_id, gRNA_id) %>% duplicated() %>% which()
resampling_results_xie <- resampling_results_xie[-idxs_to_remove,]
# next, use join to find the result to each gRNA-gene pair
results_with_type <- left_join(x = select(gRNA_gene_pair, type, pair_str), y = resampling_results_xie, by = "pair_str")
# finally, apply the fdr correction
all_results_annotated <- results_with_type %>% group_by(type) %>%
  mutate(adjusted_pvalue = ifelse(type == 'cis', p.adjust(p_value, 'fdr'), NA),
         rejected = ifelse(is.na(adjusted_pvalue), NA, adjusted_pvalue <= 0.1)) %>% ungroup()
write_fst(x = all_results_annotated, path = paste0(results_dir, "/all_results_annotated.fst"))
