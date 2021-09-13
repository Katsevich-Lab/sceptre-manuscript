code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

# load data
bulk_region_names <- paste0(processed_dir, "/bulk_region_names.rds") %>% readRDS() %>% dplyr::select(gRNA_id = region, enh_name = region_name)
gRNA_gene_pair <- read.fst(paste0(processed_dir, '/gRNA_gene_pairs.fst')) %>% mutate(pair_str = paste0(gRNA_id, "+", gene_id)) # 0 duplicated rows
gRNA_gene_pair <- dplyr::arrange(gRNA_gene_pair, pair_str)
odm <- ondisc::read_odm(metadata_fp = paste0(processed_dir, "/odm/metadata_final.rds"),
                 odm_fp = paste0(processed_dir, "/odm/gene_expression_matrix.odm"))
gene_df <- data.frame(gene_id = ondisc::get_feature_ids(odm), gene_name = ondisc::get_feature_names(odm))

# function
annotate_results <- function(results, no_bulk) {
  # append gene_names
  x <- dplyr::left_join(x = results, y = gene_df, by = "gene_id") %>%
    dplyr::left_join(., y = bulk_region_names, by = "gRNA_id") %>%
    dplyr::mutate(., gene_name = factor(gene_name), enh_name = factor(enh_name)) %>%
    mutate(pair_str = paste0(gRNA_id, "+", gene_id))
  if (no_bulk) {
    gRNA_gene_pair_sub <- gRNA_gene_pair %>% dplyr::filter(type != "bulk_validation")
    x <- dplyr::left_join(x = x, y = dplyr::select(gRNA_gene_pair_sub, pair_str, type),
                          by = "pair_str")
  } else {
    x <- dplyr::arrange(x, pair_str)
    all(x$pair_str == gRNA_gene_pair$pair_str)
    x <- x %>% dplyr::mutate(type = gRNA_gene_pair$type)
  }
  out <- x %>% dplyr::group_by(type) %>%
    mutate(adjusted_pvalue = ifelse(type == 'cis', p.adjust(p_value, 'fdr'), NA),
           rejected = ifelse(is.na(adjusted_pvalue), NA, adjusted_pvalue <= 0.1)) %>% ungroup()
  return(out)
}

# annotate three sets of results: sceptre, inb, and monocle nb
p_vals_sceptre <- paste0(results_dir, "/all_results.fst") %>% read.fst()
p_vals_inb <- paste0(results_dir_negative_binomial, "/all_results.fst") %>% read.fst()
p_vals_monocle <- paste0(results_dir_negative_binomial, "/monocle_nb_results.rds") %>% readRDS() %>% rename("p_value" = "p_val")

# annotate the results
sceptre_results_annotated <- annotate_results(p_vals_sceptre, FALSE)
inb_results_annotated <- annotate_results(p_vals_inb, FALSE)
monocle_results_annotated <- annotate_results(p_vals_monocle, TRUE)

# save the annotated results
write_fst(sceptre_results_annotated, paste0(results_dir, "/all_results_annotated.fst"))
write_fst(inb_results_annotated, paste0(results_dir_negative_binomial, "/all_results_annotated.fst"))
write_fst(monocle_results_annotated, paste0(results_dir_negative_binomial, "/monocle_results_annotated.fst"))
