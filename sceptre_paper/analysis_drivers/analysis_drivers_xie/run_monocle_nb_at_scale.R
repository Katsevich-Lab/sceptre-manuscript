small_example <- TRUE

code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")

library(monocle)
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/gasp_custom_functs.R"))
analysis_ready_dir <- paste0(offsite_dir, "data/xie/analysis_ready")

cds <- readRDS(paste0(analysis_ready_dir, "/monocole_obj.rds"))
pairs <- fst::read_fst(paste0(processed_dir, "/gRNA_gene_pairs.fst"))
if (small_example) pairs <- pairs[1:100,]

# set formula
reduced_formula <- "~ log_n_umis + batch + log_n_gRNA_umis"

# Loop over all gRNAs and genes; subset by gene, and regress on the corresponding gRNA using Molly's function
# res <- future.apply::future_lapply(X = seq(1, nrow(pairs)), FUN = function(i) {
res <- purrr::map_dfr(.x = seq(1, nrow(pairs)), .f = function(i) {
  if (i %% 50 == 0) print (i)
  gRNA <- as.character(pairs$gRNA_id[i])
  gene <- as.character(pairs$gene_id[i])
  # subset cds
  cds_sub <- cds[gene,]
  # set full formula
  full_formula <- paste0(reduced_formula, " + `", gRNA, "`")
  # regress
  de <- myDifferentialGeneTest(cds = cds_sub,
                                fullModelFormulaStr = full_formula,
                                reducedModelFormulaStr = reduced_formula,
                                relative_expr = TRUE,
                                cores = 1,
                                verbose = TRUE)
  out <- data.frame(gRNA_id = gRNA, gene_id = gene, p_val = de$pval, q_val = de$qval, beta = de$beta, intercept = de$intercept)
  return(out)
})

saveRDS(res, paste0(results_dir_negative_binomial, "/monocle_nb_results.rds"))
