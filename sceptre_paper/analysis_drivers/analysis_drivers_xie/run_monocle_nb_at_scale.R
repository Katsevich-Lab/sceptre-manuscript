small_example <- TRUE

code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")

library(monocle)
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/gasp_custom_functs.R"))
analysis_ready_dir <- paste0(offsite_dir, "data/xie/analysis_ready")

cds <- readRDS(paste0(analysis_ready_dir, "/monocole_obj.rds"))
pairs <- fst::read_fst(paste0(processed_dir, "/gRNA_gene_pairs.fst")) %>%
  dplyr::filter(type != "bulk_validation")
gRNAs <- pairs$gRNA_id %>% as.character() %>% unique()
if (small_example) gRNAs <- gRNAs[1:5]

# set formula
reduced_formula <- "~ log_n_umis + batch + log_n_gRNA_umis"

res <- purrr::map_dfr(.x = gRNAs, function(gRNA) {
  print(gRNA)
  # filter pairs according to gRNA; get curr genes
  curr_genes <- dplyr::filter(pairs, gRNA_id == !!gRNA) %>% dplyr::pull(gene_id) %>% as.character()
  if (small_example) curr_genes <- curr_genes[1:5]
  # subset cds according to curr_genes
  cds_sub <- cds[curr_genes,]
  # set full formula
  full_formula <- paste0(reduced_formula, " + `", gRNA, "`")
  # regress on all genes
  de <- myDifferentialGeneTest(cds = cds_sub,
                               fullModelFormulaStr = full_formula,
                               reducedModelFormulaStr = reduced_formula,
                               relative_expr = TRUE,
                               cores = 1,
                               verbose = TRUE)
  out <- data.frame(gRNA_id = gRNA, gene_id = de$id, p_val = de$pval, q_val = de$qval, beta = de$beta, intercept = de$intercept)
  return(out)  
})

saveRDS(res, paste0(results_dir_negative_binomial, "/monocle_nb_results.rds"))
