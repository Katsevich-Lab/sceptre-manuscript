code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

suppressPackageStartupMessages(library(R.matlab))
library(fst)
library(dplyr)
library(stringr)
library(purrr)
library(reshape2)
require(readr)

gRNA.gene.pair = read.fst(paste0(processed_dir, '/gRNA_gene_pairs.fst'))
gRNA_id = sort(as.character(unique(gRNA.gene.pair$gRNA_id)))
gRNA.fname = str_replace(gRNA_id, ':', '-')

for (i in 1:length(gRNA.fname)) {
  dest = paste0(raw_data_dir, '/', gRNA.fname[i], '-down_log-pval.mat')
  download.file(url = paste0('https://github.com/russellxie/Global-analysis-K562-enhancers/blob/master/Notebooks/Data/Hypergeometric_pvals/',
                             gRNA.fname[i], '-down_log-pval.mat?raw=true'),
                destfile = dest)
}

all_sequenced_genes_id <- rhdf5::h5read(file = paste0(raw_data_dir, "/GSM3722727_K562-dCas9-KRAB_5K-sgRNAs_Batch-4_1_filtered_gene_bc_matrices_h5.h5"),
                                        name = "/refgenome_hg38_CROP-Guide-MS2-2.1.0/genes")
xie_pfiles = setNames(paste0(raw_data_dir, '/', gRNA.fname, '-down_log-pval.mat'), gRNA_id)
extract_p_vals = function(p_mat) {
  p_vals <- exp(p_mat$matrix[1, ])
  names(p_vals) <- all_sequenced_genes_id  #actually this should be provided in download_data_2.R
  return(p_vals)
}
xie_pfiles_r <- xie_pfiles %>% map(readMat) %>% map(extract_p_vals)
saveRDS(xie_pfiles_r, file = paste0(processed_dir, '/raw_p_val_xie.rds')) # without selection of gRNA-gene pairs in gRNA.gene.pair

xie_pval_mat <- as.data.frame(xie_pfiles_r)
colnames(xie_pval_mat) <- gRNA_id

original_results_xie = do.call('rbind', lapply(1:nrow(gRNA.gene.pair), function(i) {
  data.frame(gRNA_id = gRNA.gene.pair$gRNA_id[i],
             gene_id = gRNA.gene.pair$gene_id[i],
             raw_p_val = xie_pval_mat[as.character(gRNA.gene.pair$gene_id[i]),
                                      as.character(gRNA.gene.pair$gRNA_id[i])],
             site_type = gRNA.gene.pair$type[i])
}))
rownames(original_results_xie) = NULL
saveRDS(original_results_xie, file = paste0(processed_dir, '/raw_pval_xie.rds'))


dest <- paste0(raw_data_dir, "/plot_annotation.txt")
download.file(url = "https://github.com/russellxie/Global-analysis-K562-enhancers/blob/master/Notebooks/Data/Annotations/plot_annotation.txt?raw=true", destfile = dest)
dest <- paste0(raw_data_dir, "/plotted_genes.csv")
download.file(url = "https://github.com/russellxie/Global-analysis-K562-enhancers/blob/master/Notebooks/Data/Annotations/plotted_genes.csv?raw=true", destfile = dest)

dest <- paste0(raw_data_dir, "/Perct_99.9_combined_cutoff.down_genes.mat")
download.file(url = "https://github.com/russellxie/Global-analysis-K562-enhancers/blob/master/Notebooks/Data/Gene_cutoff/Perct_99.9_combined_cutoff.down_genes.mat?raw=true", destfile = dest)


plot_annotation <- read_tsv(paste0(raw_data_dir, "/plot_annotation.txt"), col_names = c("idx", "gene_name", "chr", "pos", "strand", "color_idx", "chr_idx"))
plot_geneidx_df = read.table(paste0(raw_data_dir, '/plotted_genes.csv'), sep = '\t')
plot_geneidx_df = plot_geneidx_df[-1, 2]
p99.down = readMat(paste0(raw_data_dir, '/Perct_99.9_combined_cutoff.down_genes.mat'))$matrix[1, ]

ss.down = NULL
for(i in 1:length(gRNA.fname)){
  pval_list_down = readMat(paste0(raw_data_dir, '/', gRNA.fname[i], '-down_log-pval.mat'))$matrix[1, ]
  pval_list_down[is.infinite(pval_list_down)] = 0
  pval_list_down = - pval_list_down
  ss.down = cbind(ss.down, (pval_list_down[1:58381] - p99.down)/log(10))
  cat(i, '\t')
}
colnames(ss.down) = sort(gRNA_id)
rownames(ss.down) = all_sequenced_genes_id[1:58381]
saveRDS(ss.down, file = paste0(processed_dir, "/ss_down.rds"))

gene_id_cis = unique(gRNA.gene.pair$gene_id[gRNA.gene.pair$type == 'cis'])
idx_temp = match(gene_id_cis, all_sequenced_genes_id)

SS.down.cis = ss.down[idx_temp, ]

m.gene.id = match(gRNA.gene.pair$gene_id[gRNA.gene.pair$type == 'cis'], rownames(SS.down.cis))
m.gRNA.id = match(gRNA.gene.pair$gRNA_id[gRNA.gene.pair$type == 'cis'], colnames(SS.down.cis))
ss_xie_cis = data.frame(gene_id = gRNA.gene.pair$gene_id[gRNA.gene.pair$type == 'cis'],
                        gRNA_id = gRNA.gene.pair$gRNA_id[gRNA.gene.pair$type == 'cis'],
                        ss.down = SS.down.cis[cbind(m.gene.id, m.gRNA.id)])

ss_thres = sort(ss_xie_cis$ss.down, decreasing = T)[135]
ss_xie_cis$reject.down = ss_xie_cis$ss.down >= ss_thres
saveRDS(ss_xie_cis, file = paste0(processed_dir, '/ss_xie_cis.rds'))
