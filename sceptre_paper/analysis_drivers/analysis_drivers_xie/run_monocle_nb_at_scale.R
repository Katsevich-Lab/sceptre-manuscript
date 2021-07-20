args <- commandArgs(trailingOnly = TRUE)
library(monocle)
code_dir <- if (is.na(args[1])) "/Users/timbarry/research_code/sceptre-manuscript" else args[1]
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/gasp_custom_functs.R"))
cds <- readRDS(paste0(processed_dir, "/monocole_obj.rds"))
pairs <- fst::read_fst(paste0(processed_dir, "/gRNA_gene_pairs.fst"))
pairs <- pairs[1:100,]

# set formula
reduced_formula <- "~ log_n_umis + batch + log_n_gRNA_umis"

# set up parallel
# future::plan(future::multisession())

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