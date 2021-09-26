code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))

library(readxl)
library(dplyr)
library(katsevich2020)
library(fst)
library(readr)
library(plyranges)
library(GenomicRanges)
library(tidyr)

# create cis file for resampling results for Xie data
gRNA.gene.pair = read.fst(paste0(processed_dir, '/gRNA_gene_pairs.fst'))
all_results_annotated = read.fst(paste0(results_dir, "/all_results_annotated.fst")) %>% as_tibble()
resampling_results_xie_cis = all_results_annotated %>% filter(type == 'cis')

gene.mart = readRDS(paste0(processed_dir, '/gene_mart.rds'))
gRNA.mart = readRDS(paste0(processed_dir, '/gRNA_mart.rds'))

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

# Define a couple directories
results_dir_enrichment <- paste0(offsite_dir, "/results/xie/enrichment")
functional_data_dir <- paste0(offsite_dir, "/data/functional/")
if (!dir.exists(results_dir_enrichment)) dir.create(results_dir_enrichment)
if (!dir.exists(functional_data_dir)) dir.create(functional_data_dir)

# Read in the association results
resampling_results = resampling_results_xie_cis

ss_xie_cis = readRDS(file = paste0(processed_dir, '/ss_xie_cis.rds'))
original_results <- ss_xie_cis %>% select('gene_id', 'gRNA_id', 'ss.down', 'reject.down') %>% dplyr::rename(rejected = reject.down)

# Add annotations to Virtual FACS
original_results <- left_join(original_results, resampling_results[, c('gene_id', 'gRNA_id', 'gene_short_name', 'chr', 'strand', 'target_gene.start',
                                                   'target_gene.stop', 'TSS', 'target_site.start', 'target_site.stop',
                                                   'target_site.mid')], 
          by = c('gene_id' = 'gene_id', 'gRNA_id' = 'gRNA_id'))

monocle_results_xie = read.fst(paste0(results_dir_negative_binomial, "/monocle_results_annotated.fst"))
nb_results_xie = read.fst(paste0(results_dir_negative_binomial, "/all_results_annotated.fst"))

monocle_results_xie_cis = monocle_results_xie %>% filter(type == 'cis')
nb_results_xie_cis = nb_results_xie %>% filter(type == 'cis')

monocle_results = cbind(monocle_results_xie_cis, 
                                resampling_results[match(resampling_results$pair_str, monocle_results_xie_cis$pair_str),
                                                   c('gene_short_name', 'chr', 'strand', 'target_gene.start',
                                                     'target_gene.stop', 'TSS', 'target_site.start', 'target_site.stop',
                                                     'target_site.mid')])
nb_results = cbind(nb_results_xie_cis, 
                           resampling_results[match(resampling_results$pair_str, nb_results_xie_cis$pair_str),
                                              c('gene_short_name', 'chr', 'strand', 'target_gene.start',
                                                'target_gene.stop', 'TSS', 'target_site.start', 'target_site.stop',
                                                'target_site.mid')])

# ChIP-seq enrichment analysis

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
    summarise(rejected_vf = any(rejected)) %>%
    ungroup() %>%
    inner_join(resampling_results %>%
                 #filter(site_type == "cis") %>%
                 select(gRNA_id, rejected) %>%
                 group_by(gRNA_id) %>%
                 summarise(rejected_sceptre = any(rejected)),
               by = "gRNA_id") %>%
    inner_join(monocle_results %>% 
               select(gRNA_id, rejected) %>% 
                 group_by(gRNA_id) %>%
                 summarise(rejected_monocle = any(rejected)), 
               by = 'gRNA_id') %>%
    inner_join(nb_results %>%
                 select(gRNA_id, rejected) %>%
                 group_by(gRNA_id) %>%
                 summarise(rejected_nb = any(rejected)),
               by = 'gRNA_id') %>%
    unique() %>%
    mutate(rejected_sceptre_unique_nb = rejected_sceptre & !rejected_nb,
           rejected_sceptre_unique_monocle = rejected_sceptre & !rejected_monocle,
           rejected_sceptre_unique_vf = rejected_sceptre & ! rejected_vf,
           rejected_nb_unique = rejected_nb & !rejected_sceptre, 
           rejected_monocle_unique = rejected_monocle & !rejected_sceptre, 
           rejected_vf_unique = rejected_vf & !rejected_sceptre)

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
    rejected_sceptre = df_cand_enhancers$rejected_sceptre,
    rejected_vf = df_cand_enhancers$rejected_vf,
    rejected_monocle = df_cand_enhancers$rejected_monocle,
    rejected_nb = df_cand_enhancers$rejected_nb,
    rejected_sceptre_unique_vf = df_cand_enhancers$rejected_sceptre_unique_vf,
    rejected_sceptre_unique_monocle = df_cand_enhancers$rejected_sceptre_unique_monocle, 
    rejected_sceptre_unique_nb = df_cand_enhancers$rejected_sceptre_unique_nb,
    rejected_vf_unique = df_cand_enhancers$rejected_vf_unique,
    rejected_monocle_unique = df_cand_enhancers$rejected_monocle_unique, 
    rejected_nb_unique = df_cand_enhancers$rejected_nb_unique, 
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
    summarise(rejected_sceptre = mean(rejected_sceptre),
              rejected_vf = mean(rejected_vf),
              rejected_monocle = mean(rejected_monocle),
              rejected_nb = mean(rejected_nb),
              rejected_sceptre_unique_vf = mean(rejected_sceptre_unique_vf),
              rejected_sceptre_unique_monocle = mean(rejected_sceptre_unique_monocle), 
              rejected_sceptre_unique_nb = mean(rejected_sceptre_unique_nb),
              rejected_vf_unique = mean(rejected_vf_unique),
              rejected_monocle_unique = mean(rejected_monocle_unique), 
              rejected_nb_unique = mean(rejected_nb_unique)) %>%
    as_tibble()

  write_tsv(paired_fractions, sprintf("%s/TF_paired_enhancer_fractions.tsv", results_dir_enrichment))

  # compute odds ratios for old and new methods
  vf_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_vf, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  sceptre_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_sceptre, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  monocle_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_monocle, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  nb_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_nb, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  
  sceptre_unique_vf_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_sceptre_unique_vf, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  sceptre_unique_monocle_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_sceptre_unique_monocle, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  sceptre_unique_nb_enrichments = sapply(important_TFs,  function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_sceptre_unique_nb, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })

vf_unique_enrichment = sapply(important_TFs,  function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_vf_unique, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  monocle_unique_enrichment = sapply(important_TFs,  function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_monocle_unique, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })
  nb_unique_enrichment = sapply(important_TFs,  function(TF){
    enrichment = gr_cand_overlaps %>%
      filter(TF == !!TF, quintile %in% c(0, 2)) %>%
      as_tibble() %>%
      select(rejected_nb_unique, quintile) %>%
      table() %>%
      fisher.test()
    c(enrichment$estimate, enrichment$conf.int)
  })



  TF_enrichments = tibble(TF = important_TFs, vf_enrichments = vf_enrichments[1, ], sceptre_enrichments = sceptre_enrichments[1, ],
                          monocle_enrichments = monocle_enrichments[1, ], nb_enrichments = nb_enrichments[1, ],
                          sceptre_unique_vf_enrichments = sceptre_unique_vf_enrichments[1, ], 
                          sceptre_unique_monocle_enrichments = sceptre_unique_monocle_enrichments[1, ], 
                          sceptre_unique_nb_enrichments = sceptre_unique_nb_enrichments[1, ], 
                          vf_unique_enrichments = vf_unique_enrichment[1, ], monocle_unique_enrichments = monocle_unique_enrichment[1, ],
                          nb_unique_enrichments = nb_unique_enrichment[1, ]
                          ) %>%
    gather(method, enrichment, -TF) %>%
    mutate(method = factor(method,
                           levels = c("vf_enrichments", "sceptre_enrichments", "monocle_enrichments", "nb_enrichments",
                                      "sceptre_unique_vf_enrichments", "sceptre_unique_monocle_enrichments", 
                                      "sceptre_unique_nb_enrichments", "vf_unique_enrichments", "monocle_unique_enrichments", 
                                      "nb_unique_enrichments"),
                           labels = c("Virtual FACS", "SCEPTRE", "Monocle NB", "Improved NB", 
                                      "SCEPTRE unique Virtual FACS", "SCEPTRE unique Monocle NB", "SCEPTRE unique Improved NB", 
                                      "Virtual FACS unique", "Monocle NB unique", "Improved NB unique")))
  TF_enrichments$up_bound = c(vf_enrichments[3, ], sceptre_enrichments[3, ], monocle_enrichments[3, ], nb_enrichments[3, ],
                              sceptre_unique_vf_enrichments[3, ], sceptre_unique_monocle_enrichments[3, ], sceptre_unique_nb_enrichments[3, ],
                              vf_unique_enrichment[3, ], monocle_unique_enrichment[3, ], nb_unique_enrichment[3, ])
  TF_enrichments$lower_bound = c(vf_enrichments[2, ], sceptre_enrichments[2, ], monocle_enrichments[2, ], nb_enrichments[2, ],
                                 sceptre_unique_vf_enrichments[2, ], sceptre_unique_monocle_enrichments[2, ], 
                                 sceptre_unique_nb_enrichments[2, ],
                                 vf_unique_enrichment[2, ], monocle_unique_enrichment[2, ], nb_unique_enrichment[2, ])
  write_tsv(TF_enrichments, sprintf("%s/TF_enrichments.tsv", results_dir_enrichment))



########  HI-C enrichment analysis
  domains = read_tsv(sprintf("%sHIC/GSE63525_K562_Arrowhead_domainlist.txt", functional_data_dir),
                     col_types = "ciiciicddddd") %>%
    mutate(chr1 = sprintf("chr%s", chr1), chr2 = sprintf("chr%s", chr2))

  all_pairs = original_results %>%
    #filter(site_type == "cis") %>%
    select(chr, gene_id, target_gene.start, target_gene.stop, TSS,
           gRNA_id, target_site.start, target_site.stop, rejected) %>%
    dplyr::rename(rejected_vf = rejected) %>%
    left_join(resampling_results %>%
                #filter(site_type == "cis") %>%
                select(gene_id,  gRNA_id, rejected) %>%
                dplyr::rename(rejected_sceptre = rejected),
              by = c("gene_id", "gRNA_id")) %>%
    left_join(monocle_results %>%
                select(gene_id, gRNA_id, rejected) %>%
                dplyr::rename(rejected_monocle = rejected), 
              by = c("gene_id", "gRNA_id")) %>%
    left_join(nb_results %>%
                select(gene_id, gRNA_id, rejected) %>%
                dplyr::rename(rejected_nb = rejected), 
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

  rejected_pairs = all_pairs %>% filter(rejected_vf | rejected_sceptre | rejected_monocle | rejected_nb)
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
                      col_names = "normalization", col_types = "c") %>%
      pull() %>% as.numeric()

    rejected_pairs_chr = rejected_pairs %>%
      filter(chr == !!chr, !is.na(TAD_left)) %>%
      mutate(enhancer = 0.5*(target_site.start + target_site.stop)) %>%
      select(enhancer, TSS, TAD_left, TAD_right, gene_id, gRNA_id, rejected_vf, rejected_sceptre, rejected_monocle, rejected_nb) %>%
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

