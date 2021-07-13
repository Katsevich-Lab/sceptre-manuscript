args <- c("/Users/timbarry/research_code/sceptre-manuscript", "/Users/timbarry/research_offsite/sceptre")
code_dir <- args[1]
offsite_dir <- args[2]

library(fst)
library(magrittr)
library(sceptre)
library(monocle)
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/sceptre_function_args.R"))
gene_names <- readRDS(paste0(processed_dir, "/ordered_genes.RDS"))

gene_idx <- 1000
gene_id <- ordered_gene_ids
gene_name <- gene_names[gene_idx]
gene_exp <- cell_gene_expression_matrix[,1000]
grna_indicators <- as.integer(read_fst(gRNA_indicator_matrix_fp, "chr7:121006861-121007261")[,1])
confounders <- covariate_matrix
cell_names <- paste0("cell", seq(1, nrow(cell_gene_expression_matrix)))
exp_mat <- cell_gene_expression_matrix[,seq(1, ncol(cell_gene_expression_matrix))]

diff_test_monocle = function(gene_exp, grna_indicators, confounders, Size_Factors, disp_coeffs, 
                             gene_id, gene_name, cell_names) {
  num_cells <- length(gene_exp)
  
  exp_mat_df <- matrix(gene_exp, 1, num_cells, 
                      dimnames = list(gene_id, cell_names))
  
  pd_df <- confounders
  confounder_names <- colnames(pd_df)
  pd_df$grna <- grna_indicators
  rownames(pd_df) <- cell_names
  
  # fd_df <- data.frame(id = gene_id, gene_short_name = gene_short_name)
  fd_df <- data.frame(id = gene_id, gene_short_name = gene_name)
  rownames(fd_df) <- gene_id
  
  pd <- new("AnnotatedDataFrame", data = pd_df)
  fd <- new("AnnotatedDataFrame", data = fd_df)
  cds <- newCellDataSet(exp_mat_df,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size(),
                        lowerDetectionLimit = 0.5)
  if (FALSE) { # if known, use this code
    sizeFactors(cds) = Size_Factors
    cds@dispFitInfo$blind = c()
    
    disp_func = function(q){coefs[1] + coefs[2] / q}
    attr(disp_func, "coefficients") = data.frame(asymptDisp = disp_coeffs[1], extraPois = disp_coeffs[2])
    environment(disp_func)$coefs = disp_coeffs
    
    cds@dispFitInfo$blind$disp_func = disp_func
  }
  # else, compute from scratch
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  
  full_model_formula = sprintf("~%s", paste0(c("grna", confounder_names), collapse = "+"))
  reduced_model_formula = sprintf("~%s", paste0(confounder_names, collapse = "+"))
  diff_test_res <- differentialGeneTest(cds,
                                        fullModelFormulaStr = full_model_formula,
                                        reducedModelFormulaStr = reduced_model_formula)
  
  return(diff_test_res$pval)
}
