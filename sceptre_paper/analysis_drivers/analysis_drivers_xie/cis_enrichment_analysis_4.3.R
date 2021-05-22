args <- commandArgs(trailingOnly = TRUE) 
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE-manuscript/SCEPTRE/" else args[1] 
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

require(readxl)
library(dplyr)
require(katsevich2020)
library(fst)
library(readr)
require(plyranges)
require(GenomicRanges)
library(tidyr)

# Define a couple directories
results_dir_enrichment <- paste0(offsite_dir, "/results/xie/enrichment")
functional_data_dir <- paste0(offsite_dir, "/data/functional/")

# Read in the association results 
resampling_results_xie_cis <- paste0(results_dir, "/resampling_results_xie_cis.fst") %>% read.fst()
resampling_results = resampling_results_xie_cis
ss_xie_cis = readRDS(file = paste0(results_dir, '/ss_xie_cis.rds'))
original_results <- ss_xie_cis %>% select('gene_id', 'gRNA_id', 'ss.down', 'reject.down') %>% dplyr::rename(rejected = reject.down)
original_results = cbind(original_results, resampling_results[, c('gene_short_name', 'chr', 'strand', 'target_gene.start', 
                                                                  'target_gene.stop', 'TSS', 'target_site.start', 'target_site.stop',
                                                                  'target_site.mid')])

# ChIP-seq enrichment analysis
if (!file.exists(paste0(results_dir_enrichment, "/TF_paired_enhancer_fractions.tsv")) |
    !file.exists(paste0(results_dir_enrichment, "/TF_enrichments.tsv"))) {
  # read chipseq data
  important_TFs = c("H3K27ac", "EP300", "BRD4", "GATA2", "TAL1", "TBL1XR1", "DPF2", "RNF2")
  num_TFs = length(important_TFs)
  chipseq_data = vector("list", num_TFs)
  names(chipseq_data) = important_TFs
  for (TF in important_TFs) {
    filename = sprintf("%sChIP-seq/%s.bed", functional_data_dir, TF)
    data = suppressMessages(read_tsv(filename, col_names = FALSE))
    if (ncol(data) %in% c(8,10)){
      data = data[,c(1,2,3,7)]
    } else if(ncol(data) %in% c(5,6)){
      data = data[,c(1,2,3,5)]
    } else {
      cat(sprintf("Skipping %s because score not provided.\n", TF))
      next
    }
    names(data) = c("chrom", "chromStart", "chromEnd", "signalValue")
    data$TF = TF
    chipseq_data[[TF]] = data
  }
  chipseq_data = do.call("rbind", chipseq_data)
  
  # extract which enhancers are paired to genes in original and new analyses
  df_cand_enhancers = original_results %>%
    #filter(site_type == "cis") %>%
    select(chr, gRNA_id, target_site.start, target_site.stop, rejected) %>%
    group_by(chr, gRNA_id, target_site.start, target_site.stop) %>%
    summarise(rejected_old = any(rejected)) %>%
    ungroup() %>%
    inner_join(resampling_results %>%
                 #filter(site_type == "cis") %>%
                 select(gRNA_id, rejected) %>%
                 group_by(gRNA_id) %>%
                 summarise(rejected_new = any(rejected)),
               by = "gRNA_id") %>%
    unique() %>%
    mutate(rejected_new_unique = rejected_new & !rejected_old, 
           rejected_old_unique = rejected_old & !rejected_new)
  
  # GRanges object for chipseq data
  gr_chipseq <- GRanges(
    seqnames = chipseq_data$chrom,
    ranges = IRanges(start = chipseq_data$chromStart, end = chipseq_data$chromEnd),
    score = chipseq_data$signalValue,
    TF = chipseq_data$TF)
  
  # GRanges object for candidate enhancers
  gr_cand = GRanges(
    seqnames = df_cand_enhancers$chr,
    ranges = IRanges(start = df_cand_enhancers$target_site.start, end = df_cand_enhancers$target_site.stop),
    rejected_old = df_cand_enhancers$rejected_old,
    rejected_new = df_cand_enhancers$rejected_new,
    rejected_new_unique = df_cand_enhancers$rejected_new_unique,
    rejected_old_unique = df_cand_enhancers$rejected_old_unique,
    gRNA_id = df_cand_enhancers$gRNA_id)
  
  # Split chipseq data into quintiles
  gr_chipseq_quintiles = gr_chipseq %>%
    subsetByOverlaps(gr_cand) %>%
    group_by(TF) %>%
    mutate(quintile = 1+floor((2-1e-10)*percent_rank(score))) %>%
    ungroup()
  
  # Compute which quintile each candidate enhancer falls into
  assign_quintiles = function(TF) {
    join_overlap_left(gr_cand, gr_chipseq_quintiles %>%
                        filter(TF == !!TF)) %>%
      unique() %>%
      mutate(TF = !!TF,
             quintile = ifelse(is.na(quintile), 0, quintile))
  }
  TFs = chipseq_data %>% pull(TF) %>% unique()
  gr_cand_overlaps = do.call("c", lapply(TFs, assign_quintiles))
  
  # compute paired fractions in each quintile
  paired_fractions = gr_cand_overlaps %>%
    group_by(TF, quintile) %>%
    summarise(rejected_old = mean(rejected_old),
              rejected_new = mean(rejected_new),
              rejected_old_unique = mean(rejected_old_unique),
              rejected_new_unique = mean(rejected_new_unique)) %>%
    as_tibble()
  write_tsv(paired_fractions, sprintf("%s/TF_paired_enhancer_fractions.tsv", results_dir_enrichment))
  
  # compute odds ratios for old and new methods
  old_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_old, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  new_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_new, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  old_unique_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_old_unique, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  new_unique_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_new_unique, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  
  
  TF_enrichments = tibble(TF = important_TFs, old_enrichments = old_enrichments[1, ], new_enrichments = new_enrichments[1, ], 
                          old_unique_enrichments = old_unique_enrichments[1, ], new_unique_enrichments = new_unique_enrichments[1, ]) %>%
    gather(method, enrichment, -TF) %>%
    mutate(method = factor(method,
                           levels = c("old_enrichments", "new_enrichments",
                                      "old_unique_enrichments", "new_unique_enrichments"),
                           labels = c("Original", "Proposed", "Original Unique", "Proposed Unique")))
  TF_enrichments$up_bound = c(old_enrichments[3, ], new_enrichments[3, ], old_unique_enrichments[3, ], new_unique_enrichments[3, ])
  TF_enrichments$lower_bound = c(old_enrichments[2, ], new_enrichments[2, ], old_unique_enrichments[2, ], new_unique_enrichments[2, ])
  write_tsv(TF_enrichments, sprintf("%s/TF_enrichments.tsv", results_dir_enrichment))
}


########  HI-C enrichment analysis
if (!file.exists(sprintf("%s/rejected_pairs_HIC.tsv", results_dir_enrichment))) {
  domains = read_tsv(sprintf("%sHIC/GSE63525_K562_Arrowhead_domainlist.txt", functional_data_dir),
                     col_types = "ciiciicddddd") %>%
    mutate(chr1 = sprintf("chr%s", chr1), chr2 = sprintf("chr%s", chr2))
  
  all_pairs = original_results %>%
    #filter(site_type == "cis") %>%
    select(chr, gene_id, target_gene.start, target_gene.stop, TSS,
           gRNA_id, target_site.start, target_site.stop, rejected) %>%
    dplyr::rename(rejected_old = rejected) %>%
    left_join(resampling_results %>%
                #filter(site_type == "cis") %>%
                select(gene_id,  gRNA_id, rejected) %>%
                dplyr::rename(rejected_new = rejected),
              by = c("gene_id", "gRNA_id"))
  
  all_enhancers = all_pairs %>% select(gRNA_id, chr, target_site.start, target_site.stop) %>% unique()
  all_genes = all_pairs %>% select(gene_id, chr, target_gene.start, target_gene.stop, TSS) %>% unique()
  
  gr_enhancers = GRanges(
    seqnames = all_enhancers$chr,
    ranges = IRanges(start = all_enhancers$target_site.start, end = all_enhancers$target_site.stop),
    gRNA_id = all_enhancers$gRNA_id)
  
  gr_genes = GRanges(
    seqnames = all_genes$chr,
    ranges = IRanges(start = all_genes$target_gene.start, end = all_genes$target_gene.stop),
    gene_id = all_genes$gene_id)
  
  gr_domains <- GRanges(
    seqnames = domains$chr1,
    ranges = IRanges(start = domains$x1, end = domains$x2),
    domain_id = 1:nrow(domains))
  
  gr_genes = gr_genes %>% join_overlap_left(gr_domains)
  gr_enhancers = gr_enhancers %>% join_overlap_left(gr_domains)
  
  rejected_pairs = all_pairs %>% filter(rejected_old | rejected_new)
  num_rejected_pairs = nrow(rejected_pairs)
  TAD_left = integer(num_rejected_pairs)
  TAD_right = integer(num_rejected_pairs)
  for (pair_idx in 1:num_rejected_pairs) {
    print(pair_idx)
    overlapping_domains = intersect(gr_enhancers %>%
                                      filter(gRNA_id == rejected_pairs$gRNA_id[pair_idx]) %>%
                                      as_tibble() %>%
                                      pull(domain_id),
                                    gr_genes %>%
                                      filter(gene_id == rejected_pairs$gene_id[pair_idx]) %>%
                                      as_tibble() %>%
                                      pull(domain_id))
    merged_domain = gr_domains %>% filter(domain_id %in% overlapping_domains) %>% GenomicRanges::reduce() %>% as_tibble()
    if (nrow(merged_domain) > 0) {
      TAD_left[pair_idx] = merged_domain$start
      TAD_right[pair_idx] = merged_domain$end
    } else {
      TAD_left[pair_idx] = NA
      TAD_right[pair_idx] = NA
    }
  }
  rejected_pairs$TAD_left = TAD_left
  rejected_pairs$TAD_right = TAD_right
  
  quality = "MAPQG0"
  resolution = 5000
  resolution_name = "5kb"
  chrs = rejected_pairs %>% filter(!is.na(TAD_left)) %>% pull(chr) %>% unique() %>% sort()
  rejected_pairs_chr_list = vector("list", length(chrs))
  names(rejected_pairs_chr_list) = chrs
  
  for (chr in chrs) {
    cat(sprintf("Working on %s...\n", chr))
    cat(sprintf("Reading HI-C data...\n"))
    observed = read_tsv(sprintf("%sHIC/GSE63525_K562_intrachromosomal_contact_matrices/K562/%s_resolution_intrachromosomal/%s/%s/%s_%s.RAWobserved",
                                functional_data_dir, resolution_name, chr, quality, chr, resolution_name),
                        col_names = c("Start1", "Start2", "count"),
                        col_types = "iid")
    
    KRnorm = read_tsv(sprintf("%s/HIC/GSE63525_K562_intrachromosomal_contact_matrices/K562/%s_resolution_intrachromosomal/%s/%s/%s_%s.KRnorm",
                              functional_data_dir,resolution_name,chr,quality,chr,resolution_name),
                      col_names = "normalization", col_types = "d") %>%
      pull()
    
    rejected_pairs_chr = rejected_pairs %>%
      filter(chr == !!chr, !is.na(TAD_left)) %>%
      mutate(enhancer = 0.5*(target_site.start + target_site.stop)) %>%
      select(enhancer, TSS, TAD_left, TAD_right, gene_id, gRNA_id, rejected_old, rejected_new) %>%
      mutate_at(c("enhancer", "TSS", "TAD_left", "TAD_right"), ~floor(./resolution)+1)
    
    observed_normalized = observed %>%
      mutate(Start1 = Start1/resolution+1,
             Start2 = Start2/resolution+1,
             score = count/(KRnorm[Start1]*KRnorm[Start2])) %>%
      select(Start1, Start2, score)
    
    num_rejected_pairs = nrow(rejected_pairs_chr)
    score_ranks = numeric(num_rejected_pairs)
    
    for (idx in 1:num_rejected_pairs) {
      cat(sprintf("Working on rejected pair %d out of %d...\n", idx, num_rejected_pairs))
      enhancer = rejected_pairs_chr$enhancer[idx]
      TSS = rejected_pairs_chr$TSS[idx]
      TAD_left = rejected_pairs_chr$TAD_left[idx]
      TAD_right = rejected_pairs_chr$TAD_right[idx]
      pair_left = min(enhancer, TSS)
      pair_right = max(enhancer, TSS)
      temp = observed_normalized %>%
        filter(Start1 >= TAD_left, Start2 <= TAD_right, Start2 - Start1 == pair_right - pair_left) %>%
        mutate(score_rank = percent_rank(score)) %>%
        filter(Start1 == pair_left, Start2 == pair_right) %>%
        pull(score_rank)
      if(length(temp) == 0){
        score_ranks[idx] = NA
      }else{
        score_ranks[idx] = temp
      }
    }
    rejected_pairs_chr$score_rank = score_ranks
    rejected_pairs_chr_list[[chr]] = rejected_pairs_chr
    rm(observed)
    gc()
  }
  
  rejected_pairs_chr = do.call("rbind", rejected_pairs_chr_list)
  
  rejected_pairs = rejected_pairs %>%
    left_join(rejected_pairs_chr %>%
                select(gene_id, gRNA_id, score_rank),
              by = c("gene_id", "gRNA_id"))
  write_tsv(rejected_pairs, sprintf("%s/rejected_pairs_HIC.tsv", results_dir_enrichment))
}
