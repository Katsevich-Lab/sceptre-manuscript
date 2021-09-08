code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
analysis_ready_dir <- paste0(offsite_dir, "data/xie/analysis_ready")
library(ondisc)

# Load the necessary data
gRNA_mat <- fst::read_fst(path = paste0(analysis_ready_dir, "/gRNA_indicator_matrix.fst"))
gene_gRNA_pairs <- fst::read_fst(paste0(processed_dir, "/gRNA_gene_pairs.fst"))
gene_ids_used <- unique(gene_gRNA_pairs$gene_id) %>% as.character()
covariate_matrix <- fst::read_fst(paste0(analysis_ready_dir, "/covariate_model_matrix.fst"))

# Obtain the expression matrix
odm_fp <- paste0(processed_dir, "/odm/gene_expression_matrix.odm")
metadata_fp <- paste0(processed_dir, "/odm/metadata_final.rds")
odm <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)
exp_mat <- odm[[,seq(1, ncol(odm))]]

# ensure the dimensions of expression matrix, covariate matrix, and gRNA matrix match
dim(exp_mat); dim(covariate_matrix); dim(gRNA_mat)

# obtain the global covariate matrix
global_covariate_matrix <- cbind(covariate_matrix, gRNA_mat)

# obtain the features data frame
feature_df <- data.frame(id = get_feature_ids(odm),
                         gene_short_name = get_feature_names(odm))

# assign column and row names to the data frames and matrices
# cell names first
cell_names <- get_cell_barcodes(odm)
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
saveRDS(object = cds, file = paste0(analysis_ready_dir, "/monocole_obj.rds"))
