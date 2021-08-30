code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

#################################################
# 1. Load gene data; create gene covariate matrix
#################################################
barcodes_gene <- readRDS(paste0(processed_dir, "/cell_barcodes_gene.rds"))

# Compute the number of UMIs in each cell
exp_mat_t <- readRDS(paste0(processed_dir, "/exp_mat_t_metadata.rds")) %>% sceptre::load_fbm()
n_umis_per_cell <- big_apply(exp_mat_t, function(X, ind) {colSums(X[,ind])}) %>% unlist()
gene_covariate_matrix <- data.frame(cell_barcode = barcodes_gene, log_n_umis = n_umis_per_cell)
gene_barcode_original_order <- gene_covariate_matrix$cell_barcode

# Compute the mean expression of each gene
gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
exp_mat <- readRDS(paste0(processed_dir, "/exp_mat_metadata.rds")) %>% load_fbm()
gene_expression_p <- big_apply(exp_mat, function(X, ind) colMeans(X[,ind] >= 1)) %>% unlist()
highly_expressed_genes <- gene_ids[which(gene_expression_p >= 0.08)]

##########################################
# Load gRNA data and gRNA covariate matrix
##########################################
gRNA_indic_matrix <- readRDS(paste0(processed_dir, "/gRNA_indicator_matrix_unordered.rds"))
gRNA_covariate_matrix <- readRDS(paste0(processed_dir, "/cell_covariate_matrix_gRNA.rds"))
rownames(gRNA_covariate_matrix) <- gRNA_covariate_matrix$cell_barcode

# Intersect the cells
cells_intersect <- intersect(x = gene_barcode_original_order, gRNA_covariate_matrix$cell_barcode)

# subset the gRNA indicator matrix and covariate matrix appropriately
gRNA_indic_matrix_sub <- as.data.frame(gRNA_indic_matrix[,cells_intersect] %>% t())
gRNA_covariate_matrix_sub <- gRNA_covariate_matrix[cells_intersect,]

# get the cell subset ordering
cell_subset <- match(cells_intersect, gene_barcode_original_order)

################################
# get the global covariate matrix
################################
global_covariate_matrix <- data.frame(log_n_umis = log(n_umis_per_cell[cell_subset]),
                                      log_n_gRNA_umis = log(gRNA_covariate_matrix_sub$cell_gRNA_umi_counts),
                                      batch = gRNA_covariate_matrix_sub$batch)

# save
# create a new directory, "analysis_ready" to store all data in analysis-ready format.
analysis_ready_dir <- paste0(offsite_dir, "data/xie/analysis_ready")
if (!dir.exists(analysis_ready_dir)) dir.create(path = analysis_ready_dir, recursive = TRUE)

fst::write_fst(x = global_covariate_matrix,
               path = paste0(analysis_ready_dir, "/covariate_model_matrix.fst"))
saveRDS(object = highly_expressed_genes, file = paste0(analysis_ready_dir, "/highly_expressed_genes.rds"))
gRNA_indic_mat <- fst::write_fst(gRNA_indic_matrix_sub, paste0(analysis_ready_dir, "/gRNA_indicator_matrix.fst"))

###############################################
# save subsetted cell-by-gene expression matrix
###############################################
gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
exp_mat <- readRDS(paste0(processed_dir, "/exp_mat_metadata.rds")) %>% load_fbm()
exp_mat_mem <- exp_mat[,seq(1, ncol(exp_mat))]
exp_mat_sub <- exp_mat_mem[cell_subset,]
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

