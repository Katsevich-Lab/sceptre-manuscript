# create monocole object using Xie data
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/research_code/sceptre-manuscript" else args[1]
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

# Load the necessary data
highly_expressed_genes <- readRDS(paste0(processed_dir, "/highly_expressed_genes.rds"))
gene_ids <- as.character(readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS")))
gene_names <- as.character(readRDS(paste0(processed_dir, "/ordered_genes.RDS")))
cell_gene_expression_matrix_info <- readRDS(paste0(processed_dir, "/exp_mat_metadata.rds"))
cell_gene_expression_matrix_info$backingfile <- paste0(processed_dir, "/expression_matrix")
cell_gene_expression_matrix <- cell_gene_expression_matrix_info %>% load_fbm
gRNA_mat <- fst::read_fst(path = paste0(processed_dir, "/gRNA_indicator_matrix.fst"))
cell_subset <- readRDS(paste0(processed_dir, "/cell_subsets.rds"))[["all_cells"]]
covariate_matrix <- read.fst(paste0(processed_dir, "/covariate_model_matrix.fst"))

# Obtain the expression matrix
gene_idxs <- match(x = highly_expressed_genes, table = gene_ids)
exp_mat <- Matrix::Matrix(data = cell_gene_expression_matrix[,gene_idxs], sparse = TRUE)
exp_mat <- Matrix::t(exp_mat)

# obtain the global covariate matrix
global_covariate_matrix <- cbind(covariate_matrix, gRNA_mat)

# obtain the features data frame
feature_df <- data.frame(id = gene_ids[gene_idxs],
                         gene_short_name = gene_names[gene_idxs])

# subset exp_mat and global_covariate_matrix according to cell_subset
global_covariate_matrix <- global_covariate_matrix[cell_subset,]
exp_mat <- exp_mat[,cell_subset]

# assign column and row names to the data frames and matrices
# cell names first
cell_names <- paste0("cell", seq(1, length(cell_subset)))
colnames(exp_mat) <- cell_names
row.names(global_covariate_matrix) <- cell_names
# feature names second
feature_names <- feature_df$id
row.names(feature_df) <- feature_names
row.names(exp_mat) <- feature_names

# create annotated data frames for the cell- and feature-specific covariates
library(monocle)
cell_specific_covariates <- new("AnnotatedDataFrame", data = global_covariate_matrix)
feature_specific_covariates <- new("AnnotatedDataFrame", data = feature_df)

# create the monocole object
cds <- newCellDataSet(cellData = exp_mat,
                      phenoData = cell_specific_covariates,
                      featureData = feature_specific_covariates,
                      expressionFamily = negbinomial.size(),
                      lowerDetectionLimit = 0.5)
rm(exp_mat); rm(feature_specific_covariates); rm(cell_specific_covariates); rm(global_covariate_matrix); rm(feature_df)
gc()

# estimate the size factors and dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# save the monocle object
saveRDS(object = cds, file = paste0(processed_dir, "/monocole_obj.rds"))
