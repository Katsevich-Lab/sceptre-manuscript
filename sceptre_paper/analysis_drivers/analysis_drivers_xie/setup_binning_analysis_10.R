args <- commandArgs(trailingOnly = TRUE)
offsite_dir <- if (is.na(args[2])) .get_config_path("LOCAL_SCEPTRE_DATA_DIR") else args[1]
processed_dir <- paste0(offsite_dir, "/data/xie/processed")

# source args
# Load the expression FBM to compute the gene expressions
exp_mat_metadata <- readRDS(paste0(processed_dir, "/exp_mat_metadata.rds"))
ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
# update backing file
exp_mat_metadata$backingfile <- paste0(processed_dir, "/expression_matrix")
# load FBM given backing file
load_fbm <- function(fbm_metadata) {
  out <- bigstatsr::FBM(nrow = fbm_metadata$nrow, ncol = fbm_metadata$ncol,
                        type = fbm_metadata$type, backingfile = fbm_metadata$backingfile,
                        is_read_only = TRUE, create_bk = FALSE)
  return(out)
}
exp_mat <- load_fbm(exp_mat_metadata)

# compute the count of each gene
gene_expressions <- bigstatsr::big_apply(X = exp_mat, a.FUN = function(x, ind) {
  colSums(x[,ind])
}, a.combine = "c")

# Finally, put the gene expression counts into a data frame
gene_expressions <- data.frame(gene_id = ordered_gene_ids, expression = gene_expressions)
saveRDS(gene_expressions, paste0(processed_dir, "/gene_expressions.rds"))
