args <- commandArgs(trailingOnly = TRUE) 
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE-manuscript/SCEPTRE/" else args[1] 
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

suppressPackageStartupMessages(library(R.matlab))
library(fst)
library(dplyr)
library(stringr)
library(purrr)
library(reshape2)

gRNA.gene.pair = read.fst(paste0(processed_dir, '/gRNA_gene_pairs.fst'))
gRNA_id = as.character(unique(gRNA.gene.pair$gRNA_id))
gRNA.fname = sort(str_replace(gRNA_id, ':', '-'))

for(i in 1:length(gRNA.fname)){
  dest = paste0(raw_data_dir, '/', gRNA.fname[i], '-down_log-pval.mat')
  download.file(url = paste0('https://github.com/russellxie/Global-analysis-K562-enhancers/blob/master/Notebooks/Data/Hypergeometric_pvals/', 
                             gRNA.fname[i], '-down_log-pval.mat?raw=true'), 
                destfile = dest)
}

xie_pfiles = setNames(paste0(raw_data_dir, '/', gRNA.fname, '-down_log-pval.mat'), gRNA_id)
extract_p_vals = function(p_mat){
  p_vals <- exp(p_mat$matrix[1, ])
  names(p_vals) <- all_sequenced_genes_ids  #actually this should be provided in download_data_2.R
  return(p_vals)
}
xie_pfiles_r <- xie_pfiles %>% map(readMat) %>% map(extract_p_vals)
saveRDS(xie_pfiles_r, file = paste0(processed_dir, '/raw_p_val_xie.rds')) # without selection of gRNA-gene pairs in gRNA.gene.pair

xie_pval_mat <-do.call('cbind', xie_pfiles_r)
colnames(xie_pval_mat) <- gRNA_id
rownames(xie_pval_mat) <- all_sequenced_genes_ids

original_results_xie = do.call('rbind', lapply(1:nrow(gRNA.gene.pair), function(i){
  data.frame(gRNA_id = gRNA.gene.pair$gRNA_id[i], 
             gene_id = gRNA.gene.pair$gene_id[i], 
             raw_p_val = xie_pval_mat[gRNA.gene.pair$gene_id[i], gRNA.gene.pair$gRNA_id[i]], 
             site_type = gRNA.gene.pair$type[i])
}))
rownames(original_results_xie) = NULL

#################################
# Calcluate Significant Score 
#################################
gene.mart = readRDS(paste0(processed_dir, '/gene_mart.rds'))
gRNA.mart = readRDS(paste0(processed_dir, '/gRNA_mart.rds'))

gRNA.gene.pair = read.fst(paste0(processed_dir, '/gRNA_gene_pairs.fst'))
gRNA_id = as.character(unique(gRNA.gene.pair$gRNA_id))
gRNA.fname = sort(str_replace(gRNA_id, ':', '-'))
plot_annotation <- read_tsv(paste0(raw_data_dir, "/plot_annotation.txt"), col_names = c("idx", "gene_name", "chr", "pos", "strand",
                                                                                      "color_idx", "chr_idx"))

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
ss.down[ss.down < 0] = 0

ordered_gene_ids <- readRDS(file = paste0(processed_dir, "/ordered_gene_ids.RDS"))
ordered_genes <- readRDS(file =  paste0(processed_dir, "/ordered_genes.RDS"))

resampling_results_xie_cis = paste0(results_dir, "/resampling_results_xie_cis.fst") %>% read.fst() %>% as_tibble()
gene_id_cis = unique(resampling_results_xie_cis$gene_id)
idx_temp = plot_annotation$idx[match(ordered_genes[match(gene_id_cis, ordered_gene_ids)], plot_annotation$gene_name)]+1
idx_temp2 = match(gene_id_cis, all_sequenced_genes_id)
sum(idx_temp != idx_temp2) #0

SS.down.cis = ss.down[idx_temp, ]
SS.up.cis = ss.up[idx_temp, ]
SS.comb.cis = ss.comb[idx_temp, ]

m.gene.id = match(resampling_results_xie_cis$gene_id, rownames(SS.down.cis))
m.gRNA.id = match(resampling_results_xie_cis$gRNA_id, colnames(SS.down.cis))
ss_xie_cis = data.frame(gene_id = resampling_results_xie_cis$gene_id, 
                        gRNA_id = resampling_results_xie_cis$gRNA_id, 
                        ss.down = SS.down.cis[cbind(m.gene.id, m.gRNA.id)], 
                        ss.up = SS.up.cis[cbind(m.gene.id, m.gRNA.id)],
                        ss.comb = SS.comb.cis[cbind(m.gene.id, m.gRNA.id)])
  
ss_xie_cis$reject.down = ss_xie_cis$ss.down > 1.9
ss_xie_cis$reject.comb = ss_xie_cis$ss.down > 1.9 | ss_xie_cis$ss.up > 2.8
