code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")

# Pre-process data
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
packs <- c("bigstatsr", "future", "furrr", "Matrix", "rhdf5", "stringr", "openxlsx", "katsevich2020", "magrittr")
for (pack in packs) suppressPackageStartupMessages(library(pack, character.only = TRUE))

###################
# Expression matrix
###################
# First, we determine the number of cells and number of genes across all the batches
h5_files <- paste0(raw_data_dir, "/", grep(pattern = '*.h5', list.files(raw_data_dir), value = TRUE))
dims_across_h5s <- sapply(h5_files, function(h5_file) {
  h5_handle <- H5Fopen(h5_file)
  dim <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/shape"
  out <- as.integer(dim)
  H5Fclose(h5_handle)
  return(out)
}) %>% t()
row.names(dims_across_h5s) <- NULL
all(dims_across_h5s[,1] == dims_across_h5s[1,1]) # Verify n genes consistent across files
n_cells_total <- sum(dims_across_h5s[,2])
n_genes_total <- dims_across_h5s[1,1]

# We will use only the protein-coding genes in this analysis; load them.
gene_df <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/Genes.xlsx"), sheet = 1)
all_protein_coding_genes <- gene_df$Gene_Symbol
rm(gene_df)
h5_file <- h5_files[1]
h5_handle <- H5Fopen(h5_file)
all_sequenced_genes <- h5_handle$"/refgenome_hg38_CROP-Guide-MS2-2.1.0/gene_names"
all_sequenced_genes_ids <-  h5_handle$"/refgenome_hg38_CROP-Guide-MS2-2.1.0/genes"
saveRDS(all_sequenced_genes_ids, file = paste0(processed_dir, '/all_sequenced_genes_id.rds'))

h5_gene_info <- data.frame(row_idx = seq(1, n_genes_total), gene_name = all_sequenced_genes, gene_id = all_sequenced_genes_ids)
h5_gene_info_protein_coding <- dplyr::filter(h5_gene_info,
                                             gene_name %in% all_protein_coding_genes)
genes_in_use <- unique(h5_gene_info_protein_coding$gene_name)
genes_in_use_ids <- h5_gene_info_protein_coding$gene_id
gene_idxs_in_use <- h5_gene_info_protein_coding$row_idx
n_genes_in_use <- length(genes_in_use_ids)
H5Fclose(h5_handle)

# Next, we create a file-backed matrix to store the transpose of the expression matrix
exp_mat_t <- FBM(nrow = n_genes_in_use, ncol = n_cells_total, type = "unsigned short", init = 0, backingfile = paste0(processed_dir, "/expression_matrix_t"), create_bk = TRUE)
exp_mat_t_metadata <- list(nrow = n_genes_in_use, ncol = n_cells_total, type = "unsigned short", backingfile = paste0(processed_dir, "/expression_matrix_t"))
saveRDS(object = exp_mat_t_metadata, file = paste0(processed_dir, "/exp_mat_t_metadata.rds"))

exp_mat <- FBM(nrow = n_cells_total, ncol = n_genes_in_use, type = "unsigned short", init = 0, backingfile = paste0(processed_dir, "/expression_matrix"), create_bk = TRUE)
exp_mat_metadata <- list(nrow = n_cells_total, ncol = n_genes_in_use, type = "unsigned short", backingfile = paste0(processed_dir, "/expression_matrix"))
saveRDS(object = exp_mat_metadata, file = paste0(processed_dir, "/exp_mat_metadata.rds"))

# We iterate through the hd5 files and write the column chunks piece-by-piece to the FBM.
n_cells_processed <- 0
ordered_cell_barcodes <- character(0)

batch_lane <- stringr::str_extract(string = h5_files, pattern = "Batch-[0-9]_[0-9]") %>%
  gsub(pattern = "Batch-", replacement = "")
for (i in seq(1, length(h5_files))) {
  print(i)
  curr_batch_lane <- batch_lane[i]
  h5_file <- h5_files[i]
  print(paste("Working on", h5_file))
  h5_handle <- H5Fopen(h5_file)
  dat <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/data"
  indices <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/indices"
  ind_ptr <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/indptr"
  shape <-  h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/shape"
  n_cols <- shape[2]
  cell_barcodes <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/barcodes"
  # do some barcode processing
  if (any(duplicated(cell_barcodes))) stop("Duplicate barcodes")
  cell_barcodes <- paste0(cell_barcodes, "_", curr_batch_lane)
  # Ensure that the maximum index pointer is less than the maximum integer value in R
  if (max(ind_ptr) >= 2147483647) stop("Index pointer exceeds R max integer value.")
  sparseM <- Matrix::sparseMatrix(i = indices, p = ind_ptr, x = dat, dims = shape)
  sparseM <- sparseM[gene_idxs_in_use,]
  # write
  exp_mat_t[seq(1, n_genes_in_use), (n_cells_processed + 1):(n_cells_processed + n_cols)] <- as.matrix(sparseM)
  exp_mat[(n_cells_processed + 1):(n_cells_processed + n_cols), seq(1, n_genes_in_use)] <- as.matrix(t(sparseM))
  H5Fclose(h5_handle)
  n_cells_processed <- n_cells_processed + n_cols
  ordered_cell_barcodes <- c(ordered_cell_barcodes, cell_barcodes)
}

saveRDS(object = ordered_cell_barcodes, file = paste0(processed_dir, "/cell_barcodes_gene.rds"))
saveRDS(object = genes_in_use, file = paste0(processed_dir, "/ordered_genes.RDS"))
saveRDS(object = genes_in_use_ids, file = paste0(processed_dir, "/ordered_gene_ids.RDS"))

###############
# gRNA UMI data
###############
rm(list=ls())
code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

# Save the bulk region names
enh_targets_df <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/enh_targets.xlsx"), sheet = 1)
bulk_region_names <- filter(enh_targets_df, gene_names %in% c("ARL15", "MYB")) %>%
  select(region, region_name = Denoted.Region.Name.used.in.the.paper, targeted_gene = gene_names) %>%
  filter(region_name %in% c("ARL15-enh", "MYB-enh-3"))
saveRDS(object = bulk_region_names, paste0(processed_dir, "/bulk_region_names.rds"))

####################################################
# Next, we create a data frame containing UMI counts
####################################################
guide_seqs <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/all_oligos.xlsx"),
                                  sheet = 1) %>% dplyr::select(hg38_enh_region = "region.pos.(hg38)",
                                                               spacer_seq = "spacer.sequence")
weird_spacers <- names(which(table(guide_seqs$spacer_seq) >= 2))
guide_seqs <- dplyr::filter(guide_seqs, !(spacer_seq %in% weird_spacers))
all(table(guide_seqs$spacer_seq) == 1); all(table(guide_seqs$hg38_enh_region) == 10)

# get the raw file names
raw_fs <- list.files(raw_data_dir)
gRNA_files <- paste0(raw_data_dir, "/", grep(pattern = "sgRNA-enrichment_5K-sgRNAs_Batch", x = raw_fs, value = TRUE))

# get the batch ID of each file
batch <- stringr::str_extract(string = gRNA_files, pattern = "Batch_[0-9]_[0-9]") %>%
  gsub(pattern = "Batch_", replacement = "")

# construct a sparse matrix of gRNA counts for each file
res <- lapply(X = seq(1L, length(gRNA_files)), FUN = function(i) {
  curr_file <- gRNA_files[i]
  curr_batch <- batch[i]
  print(paste("Working on file", curr_file))
  curr_gRNA_count_matrix <- readr::read_tsv(file = curr_file,
                                            col_names = c("cell_barcode", "total_read_count", "total_umi_count", "gRNA_spacer_seqs", "read_counts", "umi_counts"),
                                            col_types = c("cccccc")) %>% dplyr::select(cell_barcode, gRNA_spacer_seqs, umi_counts, total_umi_count)
  cell_barcodes <- curr_gRNA_count_matrix$cell_barcode
  # convert to sparse triplet matrix
  umi_count_list <- lapply(curr_gRNA_count_matrix$umi_counts, function(v) stringr::str_split(v, pattern = ";") %>% unlist() %>% as.integer())
  gRNA_spacer_seq_list <- lapply(curr_gRNA_count_matrix$gRNA_spacer_seqs, function(v) stringr::str_split(v, pattern = ";") %>% unlist())
  n_entries_per_cell <- sapply(X = umi_count_list, function(v) length(v))
  barcode_vector <- rep(x = cell_barcodes, times = n_entries_per_cell)
  gRNA_spacer_seq_vector <- unlist(gRNA_spacer_seq_list)
  umi_count_vector <- unlist(umi_count_list)
  ch_df <- data.frame(barcode = barcode_vector,
                      spacer_seq = gRNA_spacer_seq_vector,
                      umi_count = umi_count_vector)
  # remove all entries with spacer seqs not in the spacer seq list
  ch_df <- ch_df %>% dplyr::filter(gRNA_spacer_seq_vector %in% guide_seqs$spacer_seq)
  ordered_barcode_vector <- sort(unique(ch_df$barcode))
  
  spacer_seq_idx <- match(x = ch_df$spacer_seq, guide_seqs$spacer_seq)
  barcode_idx <- match(x = ch_df$barcode, table = ordered_barcode_vector)
  
  m <- Matrix::sparseMatrix(i = spacer_seq_idx,
                            j = barcode_idx,
                            x = ch_df$umi_count,
                            dims = c(length(guide_seqs$spacer_seq),
                                     length(ordered_barcode_vector)))
  rownames(m) <- guide_seqs$spacer_seq
  colnames(m) <- paste0(ordered_barcode_vector, "-1_", curr_batch)
  # obtain the gRNA UMI counts and cell count
  curr_gRNA_count_matrix <- as.data.frame(curr_gRNA_count_matrix)
  rownames(curr_gRNA_count_matrix) <- curr_gRNA_count_matrix$cell_barcode
  curr_gRNA_count_matrix_sub <- curr_gRNA_count_matrix[ordered_barcode_vector,]
  umi_counts <- as.integer(curr_gRNA_count_matrix_sub$total_umi_count)
  n_cells <- length(umi_counts)
  return(list(m = m, umi_counts = umi_counts, n_cells = n_cells))
})

combined_m <- do.call(what = cbind, args = lapply(res, function(i) i$m))
all(Matrix::colSums(combined_m) >= 1); all(Matrix::rowSums(combined_m) >= 1)

# Finally, perform the grouping operation
gRNA_matrix_raw <- combined_m
unique_regions <- unique(guide_seqs$hg38_enh_region)

thresholded_ungrouped_mat <- apply(X = gRNA_matrix_raw, MARGIN = 1, FUN = function(r) {
  r_nonzero <- r[r >= 1]
  as.integer(r > mean(r_nonzero))
}) %>% Matrix::t()
colnames(thresholded_ungrouped_mat) <- colnames(gRNA_matrix_raw)

thresholded_grouped_mat <- sapply(X = unique_regions, FUN = function(curr_region) {
  print(curr_region)
  curr_spacer_seqs <- dplyr::filter(guide_seqs, hg38_enh_region == curr_region) %>%
    dplyr::pull(spacer_seq)
  curr_m <- thresholded_ungrouped_mat[curr_spacer_seqs,]
  apply(X = curr_m, MARGIN = 2, FUN = function(col) max(col))
}) %>% t()

saveRDS(thresholded_grouped_mat, paste0(processed_dir, "/gRNA_indicator_matrix_unordered.rds"))


# compute the cell covariate matrix
cell_gRNA_umi_counts <- lapply(res, function(i) i$umi_counts) %>% do.call(what = "c", args = .)
cell_counts <- lapply(res, function(i) i$n_cells) %>% unlist()
cell_covariate_matrix <- data.frame(cell_barcode = colnames(thresholded_grouped_mat),
                                    cell_gRNA_umi_counts = cell_gRNA_umi_counts,
                                    batch = rep(batch, times = cell_counts))
saveRDS(object = cell_covariate_matrix, file = paste0(processed_dir, "/cell_covariate_matrix_gRNA.rds"))

##############
# Bulk RNA-seq
##############
library(openxlsx)
rm(list=ls())
code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

gene_df <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/Genes.xlsx"), sheet = 1)
all_protein_coding_genes <- gene_df$Gene_Symbol
bulk_info <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/bulk_rna_info.xlsx"), sheet = 3) %>% select(library_name = Library.Name, gRNA = sgRNA, region = Region, biological_duplicate = Biological.Duplicate)
bulk_info_arl15_enh <- slice(bulk_info, 1:25)
bulk_info_myb_enh3 <- slice(bulk_info, 26:49)

bulk_df <- suppressWarnings(read_tsv(file = paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.ARL15_enhancer.txt"), col_types = "ccccciiiiiiiiiiiiiiiiiiiiiiiiii")) %>% rename("PZ788" = "PZ778...10", "PZ778" = "PZ778...13")
bulk_df_arl15_enh <- filter(bulk_df, Geneid %in% all_protein_coding_genes)
bulk_df <- suppressWarnings(read_tsv(file = paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.MYB_enhancer.txt"), col_types = "ccccciiiiiiiiiiiiiiiiiiiiiiiiii"))
bulk_df_myb_enh3 <- filter(bulk_df, Geneid %in% all_protein_coding_genes)

bulk_rnaseq <- list(data = list(arl15_enh = bulk_df_arl15_enh, myb_enh3 = bulk_df_myb_enh3), info = list(arl15_enh = bulk_info_arl15_enh, myb_enh3 = bulk_info_myb_enh3))
saveRDS(object = bulk_rnaseq, file = paste0(processed_dir, "/bulk_RNAseq.rds"))

#################################################
# Xie hypergeometric p-values for ARL15 and MYB-3
#################################################
all_sequenced_genes_ids <- readRDS(file = paste0(processed_dir, '/all_sequenced_genes_id.rds'))
suppressPackageStartupMessages(library(R.matlab))
extract_p_vals_hypergeo <- function(p_vals_hypergeo) {
  p_vals <- exp(p_vals_hypergeo$matrix[1,])
  names(p_vals) <- all_sequenced_genes_ids
  return(p_vals)
}

xie_pfiles <- setNames(paste0(raw_data_dir, c("/hypergeometric_pvals_arl15_down.mat", "/hypergeometric_pvals_myb3_down.mat")), c("arl15_enh", "myb_enh3"))
xie_pfiles_r <- xie_pfiles %>% map(readMat) %>% map(extract_p_vals_hypergeo)

saveRDS(object = xie_pfiles_r, file = paste0(processed_dir, c("/xie_p_values.rds")))
