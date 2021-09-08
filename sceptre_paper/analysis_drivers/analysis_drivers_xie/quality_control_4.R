code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
library(ondisc)

##############
# 1. Gene data
##############
odm_dir <- paste0(processed_dir, "/odm")
odm_fp <- paste0(odm_dir, "/gene_expression_matrix.odm")
metadata_fp <- paste0(odm_dir, "/metadata_plus.rds")
gene_odm <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)

# Compute the mean expression of each gene
highly_expressed_genes <- get_feature_covariates(gene_odm) %>% dplyr::mutate(gene_id = get_feature_ids(gene_odm),
                                                   p_expressed = n_nonzero / ncol(gene_odm)) %>%
  dplyr::filter(p_expressed >= 0.005, protein_coding) %>% dplyr::pull(gene_id)

# Save the highly expressed genes; subset expression matrix according to these genes
saveRDS(object = highly_expressed_genes, file = paste0(processed_dir, "/highly_expressed_genes.rds"))
gene_odm <- gene_odm[highly_expressed_genes,] # 10,014 genes, 106,670 cells

##############
# 2. gRNA data
##############
gRNA_indic_matrix <- readRDS(paste0(processed_dir, "/gRNA_indicator_matrix_unordered.rds"))
gRNA_covariate_matrix <- readRDS(paste0(processed_dir, "/cell_covariate_matrix_gRNA.rds"))
rownames(gRNA_covariate_matrix) <- gRNA_covariate_matrix$cell_barcode

######################
# 3. Cell intersection
######################
# Intersect the gRNA and gene cells
cells_intersect <- intersect(x = get_cell_barcodes(gene_odm),
                             y = colnames(gRNA_indic_matrix))

# subset the gRNA indicator matrix and covariate matrix appropriately
gRNA_indic_matrix_sub <- as.data.frame(gRNA_indic_matrix[,cells_intersect] %>% t())
gRNA_covariate_matrix_sub <- gRNA_covariate_matrix[cells_intersect,]

# subset the gene matrix
gene_odm <- gene_odm[,cells_intersect]
# obtain the global cell covariate matrix
global_covariate_matrix <- data.frame(log_n_umis = get_cell_covariates(gene_odm) %>% dplyr::pull(n_umis) %>% log(),
                                      log_n_gRNA_umis = log(gRNA_covariate_matrix_sub$cell_gRNA_umi_counts),
                                      batch = gRNA_covariate_matrix_sub$batch)
save_odm(odm = gene_odm, metadata_fp = paste0(processed_dir, "/odm/metadata_final.rds"))
##################################################################################
# Create analysis-ready dir; save highly expressed genes and gRNA indicator matrix
##################################################################################
# create a new directory, "analysis_ready" to store all data in analysis-ready format.
analysis_ready_dir <- paste0(offsite_dir, "data/xie/analysis_ready")
if (!dir.exists(analysis_ready_dir)) dir.create(path = analysis_ready_dir, recursive = TRUE)

fst::write_fst(x = global_covariate_matrix,
               path = paste0(analysis_ready_dir, "/covariate_model_matrix.fst"))
fst::write_fst(gRNA_indic_matrix_sub, paste0(analysis_ready_dir, "/gRNA_indicator_matrix.fst"))
saveRDS(highly_expressed_genes, paste0(analysis_ready_dir, "/ordered_gene_ids.rds"))

######################################################################################
# save subsetted cell-by-gene expression matrix (only highly expressed genes included)
######################################################################################
exp_mat_sub <- t(as.matrix(gene_odm[[,1:ncol(gene_odm)]]))
exp_mat_sub_disk <- FBM(nrow = nrow(exp_mat_sub),
                        ncol = ncol(exp_mat_sub),
                        type = "unsigned short",
                        init = 0,
                        backingfile = paste0(analysis_ready_dir, "/expression_matrix_sub"),
                        create_bk = TRUE)
exp_mat_sub_disk[1:nrow(exp_mat_sub), 1:ncol(exp_mat_sub)] <- exp_mat_sub
exp_mat_sub_metadata <- list(nrow = nrow(exp_mat_sub),
                             ncol = ncol(exp_mat_sub),
                             type = "unsigned short",
                             backingfile = paste0(analysis_ready_dir, "/expression_matrix_sub"))
saveRDS(object = exp_mat_sub_metadata, paste0(analysis_ready_dir, "/exp_mat_sub_metadata.rds"))
