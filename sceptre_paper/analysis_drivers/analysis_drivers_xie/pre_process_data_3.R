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

genes_in_use <- readRDS(paste0(processed_dir, "/ordered_genes.RDS"))
genes_in_use_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))

# Load the gRNA identification information; we save the regions of ARL15-enh and MYB-enh-1-4.
enh_targets_df <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/enh_targets.xlsx"), sheet = 1)
bulk_region_names <- filter(enh_targets_df, gene_names %in% c("ARL15", "MYB")) %>%
  select(region, region_name = Denoted.Region.Name.used.in.the.paper, targeted_gene = gene_names) %>%
  filter(region_name %in% c("ARL15-enh", "MYB-enh-3"))
saveRDS(object = bulk_region_names, paste0(processed_dir, "/bulk_region_names.rds"))

guide_seqs <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/all_oligos.xlsx"), sheet = 1) %>% rename(hg38_enh_region = "region.pos.(hg38)")
target_regions <- guide_seqs %>% pull(hg38_enh_region) %>% unique()
names(target_regions) <- target_regions

# Next, we create a data frame containing UMI counts for all targeted putative enhancers
raw_fs <- list.files(raw_data_dir)
gRNA_files <- paste0(raw_data_dir, "/", grep(pattern = "sgRNA-enrichment_5K-sgRNAs_Batch", x = raw_fs, value = TRUE))
batch_lane <- stringr::str_extract(string = gRNA_files, pattern = "Batch_[0-9]_[0-9]") %>%
  gsub(pattern = "Batch_", replacement = "")

# Run this part in parallel
plan(multisession, workers = 10)
res <- future_map(.x = seq(1L, length(gRNA_files)), .f = function(i) {
  curr_file <- gRNA_files[i]
  curr_batch_lane <- batch_lane[i]
  print(paste("Working on file", curr_file))
  curr_gRNA_count_matrix <- read_tsv(file = curr_file, col_names = c("cell_barcode", "total_read_count", "total_umi_count", "gRNA_spacer_seqs", "read_counts", "umi_counts"), col_types = c("cccccc")) %>% select(cell_barcode, gRNA_spacer_seqs, umi_counts, total_umi_count)
  cell_barcodes <- pull(curr_gRNA_count_matrix, cell_barcode)
  # check for duplicates
  if (any(duplicated(cell_barcodes))) stop("Duplicated barcodes.")
  cell_barcodes <- paste0(cell_barcodes, "-1_", curr_batch_lane)
  
  count_matrix_list <- map(target_regions, function(region) {
    cat(paste0("Working on region ", region, ".\n"))
    region_spacer_seqs <- filter(guide_seqs, hg38_enh_region == region) %>% pull(spacer.sequence)
    curr_batch_gRNA_umi_counts <- sapply(X = 1:nrow(curr_gRNA_count_matrix), FUN = function(row_id) {
      r <- curr_gRNA_count_matrix[row_id,]
      spacers <- str_split(r$gRNA_spacer_seqs, pattern = ";") %>% unlist()
      umi_counts <- str_split(r$umi_counts, pattern = ";") %>% unlist() %>% as.integer()
      umi_locs <- match(x = region_spacer_seqs, table = spacers)
      curr_counts <- sapply(umi_locs, function(i) if (is.na(i)) 0 else umi_counts[i])
      names(curr_counts) <- region_spacer_seqs
      curr_counts
    }) %>% t() %>% Matrix(sparse = TRUE)
  })
  list(cell_barcodes = cell_barcodes, count_matrix_list = count_matrix_list, total_umis = as.integer(curr_gRNA_count_matrix$total_umi_count))
})

gRNA_count_matrix_list <- map(target_regions, function(region) {
  map(res, function(x) x$count_matrix_list[[region]]) %>% reduce(.f = rbind)
})

cell_barcodes <- map(.x = res, .f = function(x) x$cell_barcodes) %>% reduce(.f = c)
cell_counts <- sapply(X = res, FUN = function(x) length(x$cell_barcodes))
cell_gRNA_umi_counts <- map(.x = res, function(x) x$total_umis) %>% reduce(.f = c)

# We reduce each matrix in the gRNA_count_matrix_list to a single logical vector
combine_gRNAs_in_group <- function(gRNA_count_matrix) {
  apply(X = gRNA_count_matrix, MARGIN = 2, FUN = function(column) {
    v <- sum(column >= 1)
    U <- sum(column)
    column/U > 1/v
  }) %>% apply(MARGIN = 1, FUN = function(r) any(r))
}
gRNA_indic_matrix <- map_dfr(gRNA_count_matrix_list, combine_gRNAs_in_group) %>% as.data.frame()
rownames(gRNA_indic_matrix) <- cell_barcodes
saveRDS(object = gRNA_indic_matrix, file = paste0(processed_dir, "/gRNA_indicator_matrix_unordered.rds"))

# create the cell covariate gRNA matrix
cell_covariate_matrix <- data.frame(cell_barcode = cell_barcodes,
                                    cell_gRNA_umi_counts = cell_gRNA_umi_counts,
                                    batch = rep(batch_lane, times = cell_counts))
saveRDS(object = cell_covariate_matrix, file = paste0(processed_dir, "/cell_covariate_matrix_gRNA.rds"))

##############
# Bulk RNA-seq
##############
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
