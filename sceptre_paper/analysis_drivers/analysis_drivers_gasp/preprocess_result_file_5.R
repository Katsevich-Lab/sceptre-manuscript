args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
require(readxl)
require(katsevich2020)
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_gasp/file_paths_to_dirs.R"))

# Load the Gasperini results
original_results_raw <- suppressWarnings(read_tsv(paste0(raw_data_dir, "/GSE120861_all_deg_results.at_scale.txt"), col_types = "cddddddccccciiciiccc"))
old_pairs <- read_excel(path = paste0(raw_data_dir, "/Gasperini_TableS2.xlsx"), sheet = 3)

get_target_site <- function(grna_group) {
  if(!grepl("_two", grna_group)){
    target_site = grna_group
  } else{
    target_site = strsplit(grna_group, "_")[[1]][1]
  }
}

# extract target sites, high confidence subset
original_results <- original_results_raw %>%
  group_by(gRNA_group) %>%
  mutate(target_site = get_target_site(unique(gRNA_group))) %>%
  ungroup() %>%
  left_join(old_pairs %>%
              select(Target_Site, ENSG, high_confidence_subset) %>%
              mutate(quality_rank_grna = "top_two") %>%
              rename(target_site = Target_Site),
            by = c("target_site", "ENSG", "quality_rank_grna")) %>%
  rename(grna_group = gRNA_group, gene_id = ENSG, pair_id = pairs4merge, chr = target_gene.chr) %>%
  mutate(TSS = ifelse(strand == "+", target_gene.start, target_gene.stop),
         rejected = ifelse(is.na(high_confidence_subset), FALSE, high_confidence_subset)) %>%
  select(chr, pair_id, rejected,
         gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene,
         grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type,
         beta, intercept, fold_change.transcript_remaining, pvalue.raw, pvalue.empirical, pvalue.empirical.adjusted)

# fill in missing empirical p-values
null_ecdf <- ecdf(c(0,original_results %>% filter(site_type == "NTC") %>% pull(pvalue.raw)))
original_results <- original_results %>%
  mutate(pvalue.empirical = ifelse(is.na(pvalue.empirical), null_ecdf(pvalue.raw), pvalue.empirical))

# write to file
write.fst(original_results, paste0(processed_dir, "/original_results.fst"))

# Manipulate the SCEPTRE results
sceptre_res <- paste0(results_dir, "/all_results.fst") %>% read.fst()
resampling_results <- sceptre_res %>% rename(grna_group = gRNA_id) %>% left_join(original_results %>%
                                                             select(chr, pair_id,
                                                                    gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene,
                                                                    grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type),
                                                           by = c("gene_id", "grna_group")) %>%
  group_by(site_type, quality_rank_grna) %>%
  mutate(adjusted_pvalue = ifelse(site_type == "DHS" & quality_rank_grna == "top_two",
                                  p.adjust(p_value, "fdr"), NA),
         rejected = ifelse(is.na(adjusted_pvalue), FALSE, adjusted_pvalue <= 0.1)) %>%
  ungroup() %>%
  select(chr, pair_id, rejected,
         gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene,
         grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type,
         p_value, xi, omega, alpha, nu)

# save
write.fst(resampling_results, paste0(results_dir, "/resampling_results.fst"))

# Manipulate the negative binomial regression results
negbin_results <- paste0(results_dir_negative_binomial, "/all_results.fst") %>% read.fst()
likelihood_results <- negbin_results %>% rename(grna_group = gRNA_id) %>%
  left_join(original_results %>%
              select(chr, pair_id,
                     gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene,
                     grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type),
            by = c("gene_id", "grna_group")) %>%
  select(chr, pair_id,
         gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene,
         grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type,
         pvalue = p_value)

write.fst(likelihood_results, paste0(results_dir_negative_binomial, "/likelihood_results.fst"))
