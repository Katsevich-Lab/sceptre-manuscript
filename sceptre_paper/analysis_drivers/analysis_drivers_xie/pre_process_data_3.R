code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
library(ondisc)
library(rhdf5)
library(openxlsx)

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
batch_lane <- stringr::str_extract(string = h5_files, pattern = "Batch-[0-9]_[0-9]") %>%
  gsub(pattern = "Batch-", replacement = "")

# Next, initialize an ondisc matrix storing the gene expression data
if (!dir.exists(paste0(processed_dir, "/odm"))) dir.create(paste0(processed_dir, "/odm"))
h5_odm <- create_ondisc_matrix_from_h5_list(h5_list = h5_files,
                                            odm_fp = paste0(processed_dir, "/odm/gene_expression_matrix"),
                                            metadata_fp = paste0(processed_dir, "/odm/gene_expression_metadata.rds"),
                                            barcode_suffixes = batch_lane,
                                            progress = TRUE)
h5_odm <- read_odm(odm_fp = paste0(processed_dir, "/odm/gene_expression_matrix.odm"),
                   metadata_fp = paste0(processed_dir, "/odm/gene_expression_metadata.rds"))

# modify the ondisc matrix by identifying the protein-coding genes
gene_df <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/Genes.xlsx"), sheet = 1)
all_protein_coding_genes <- gene_df$Gene_Symbol
rm(gene_df)
is_protein_coding <- get_feature_names(h5_odm) %in% all_protein_coding_genes
h5_odm_w_protein_coding <- h5_odm %>% mutate_feature_covariates(protein_coding = is_protein_coding)
save_odm(h5_odm_w_protein_coding, paste0(processed_dir, "/odm/metadata_plus.rds"))
all_gene_ids <- rhdf5::h5read(file = h5_files[1], name = "/refgenome_hg38_CROP-Guide-MS2-2.1.0/genes")
saveRDS(object = all_gene_ids,
        file = paste0(processed_dir, '/all_sequenced_genes_id.rds'))

###############
# gRNA UMI data
###############
rm(list=ls())
code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

# Save the bulk region names
h5_files <- paste0(raw_data_dir, "/", grep(pattern = '*.h5', list.files(raw_data_dir), value = TRUE))
gene_names <- rhdf5::h5read(h5_files[1], "/refgenome_hg38_CROP-Guide-MS2-2.1.0/gene_names")
gene_ids <- rhdf5::h5read(h5_files[1], "/refgenome_hg38_CROP-Guide-MS2-2.1.0/genes")

enh_targets_df <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/enh_targets.xlsx"), sheet = 1)
bulk_region_names <- filter(enh_targets_df, gene_names %in% c("ARL15", "MYB")) %>%
  dplyr::select(region, region_name = Denoted.Region.Name.used.in.the.paper, targeted_gene = gene_names) %>%
  filter(region_name %in% c("ARL15-enh", "MYB-enh-3")) %>%
  dplyr::mutate(targeted_gene_id = c(gene_ids["ARL15" == gene_names],
                                     gene_ids["MYB" == gene_names]))
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
