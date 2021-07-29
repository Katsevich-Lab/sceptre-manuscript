offsite_dir <- if (is.na(args[2])) "/Volumes/tims_new_drive/research/sceptre_files/" else args[2]
manuscript_figure_dir <- paste0(code_dir, "/sceptre_paper/manuscript/figures")

# Gasperini results
original_results_gasp = paste0(offsite_dir, "/data/gasperini/processed/original_results.fst") %>% read.fst()
resampling_results_gasp = read.fst(sprintf("%s/results/gasperini/sceptre/resampling_results.fst", offsite_dir)) %>% as_tibble()
likelihood_results_gasp = read.fst(sprintf("%s/results/gasperini/negative_binomial/likelihood_results.fst", offsite_dir)) %>% as_tibble()
covariates_gasp = read.fst(sprintf("%s/data/gasperini/processed/cell_covariate_model_matrix.fst", offsite_dir)) %>% as_tibble()

# Xie results: 4 methods: original, resampling_results (sceptre), likelihood results (iNB), and monocle_nb
gRNA.gene.pair = read.fst(sprintf("%s/data/xie/processed/gRNA_gene_pairs.fst", offsite_dir)) %>% as_tibble()
original_results_xie = readRDS(sprintf("%s/data/xie/processed/raw_pval_xie.rds", offsite_dir)) %>% as_tibble()
resampling_results_xie_with_names = read.fst(sprintf("%s/results/xie/sceptre/all_results_annotated.fst", offsite_dir)) %>% as_tibble()
likelihood_results_xie = read.fst(sprintf("%s/results/xie/negative_binomial/all_results.fst", offsite_dir)) %>% as_tibble()

p_vals_bulk <- paste0(offsite_dir, "/results/xie/bulk_rna_seq/pvals_arl15_enh.rds") %>% readRDS() %>% as_tibble() %>% rename(gene_names = gene_id)
p_vals_bulk_myb3_enh3 <- paste0(offsite_dir, "/results/xie/bulk_rna_seq/pvals_myb_enh3.rds") %>% readRDS() %>% as_tibble() %>% rename(gene_names = gene_id)
covariates_xie = read.fst(sprintf("%s/data/xie/processed/covariate_model_matrix.fst", offsite_dir)) %>% as_tibble()
grna_indicator_matrix_xie = read.fst(sprintf("%s/data/xie/processed/gRNA_indicator_matrix.fst", offsite_dir)) %>% as_tibble()
ss_xie_cis = readRDS(sprintf("%s/data/xie/processed/ss_xie_cis.rds", offsite_dir)) %>% as_tibble()

resampling_results_xie_cis = resampling_results_xie_with_names %>% filter(type == 'cis')
gene.mart = readRDS(sprintf("%s/data/xie/processed/gene_mart.rds", offsite_dir))
gRNA.mart = readRDS(sprintf("%s/data/xie/processed/gRNA_mart.rds", offsite_dir))
gene.ensembl.id <- lapply(strsplit(as.character(resampling_results_xie_cis$gene_id), '[.]'), function(x){x[1]}) %>% unlist
resampling_results_xie_cis = cbind(resampling_results_xie_cis,
                                   gene.mart[match(gene.ensembl.id, gene.mart$ensembl_gene_id),
                                             c('hgnc_symbol', 'chr', 'strand', 'start_position', 'end_position')])
resampling_results_xie_cis <- resampling_results_xie_cis %>% dplyr::rename(gene_short_name = hgnc_symbol,
                                                                    target_gene.start = start_position,
                                                                    target_gene.stop = end_position)
resampling_results_xie_cis <- resampling_results_xie_cis %>% mutate(TSS = ifelse(strand > 0, target_gene.start, target_gene.stop))
resampling_results_xie_cis = cbind(resampling_results_xie_cis, gRNA.mart[match(resampling_results_xie_cis$gRNA_id, gRNA.mart$gRNA_id),
                                             c('start_position', 'end_position', 'mid_position')])
resampling_results_xie_cis <- resampling_results_xie_cis %>% dplyr::rename(target_site.start = start_position,
                                                                    target_site.stop = end_position,
                                                                    target_site.mid = mid_position)


# simulation results
simulation_results = read.fst(sprintf("%s/results/simulations/all_results.fst", offsite_dir)) %>% as_tibble()

# enrichment results Gasperini
rejected_pairs_HIC = read_tsv(sprintf("%s/results/gasperini/enrichment/rejected_pairs_HIC.tsv", offsite_dir),
                              col_types = "cciiiciilliid")
TF_enrichments = read_tsv(sprintf("%s/results/gasperini/enrichment/TF_enrichments.tsv", offsite_dir),
                          col_types = "ccd")
paired_fractions = read_tsv(sprintf("%s/results/gasperini/enrichment/TF_paired_enhancer_fractions.tsv", offsite_dir),
                            col_types = "cidddd")

# enrichment results Xie
rejected_pairs_HIC_xie <- read_tsv(sprintf("%s/results/xie/enrichment/rejected_pairs_HIC.tsv", offsite_dir),
                              col_types = "cciiiciilliid") %>% as_tibble()
TF_enrichments_xie = read_tsv(sprintf("%s/results/xie/enrichment/TF_enrichments.tsv", offsite_dir),
                              col_types = "ccddd") %>% as_tibble()
paired_fractions_xie = read_tsv(sprintf("%s/results/xie/enrichment/TF_paired_enhancer_fractions.tsv", offsite_dir),
                            col_types = "cidddd") %>% as_tibble()

# dispersion coefficients alpha_0 and alpha_1 from equation (6) in DESeq2 paper
disp_coeffs = as.numeric(readRDS(sprintf("%s/data/gasperini/processed/disp_coefficients.rds", offsite_dir)))

# dispersion table (mean expressions and dispersion estimates for each gene)
disp_table = read.fst(sprintf("%s/data/gasperini/processed/disp_table.fst", offsite_dir)) %>% as_tibble()

# additional items
ci <- 0.95
plot_colors <- c(gasperini_nb = "firebrick3", hf_nb = "hotpink2", fixed_dispersion_nb = "darkslategray4", sceptre = "royalblue4", hypergeometric = "darkslategray4", scMAGeCK = "magenta4")
subsampling_factor <- 250
