args <- commandArgs(trailingOnly = TRUE) 
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE-manuscript/SCEPTRE/" else args[1] 
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))


# Both gene ID and gRNA position are get from GRCh 38. 

#####################
# 1. gene positions 
#####################
gene.id <- readRDS(paste0(processed_dir, '/highly_expressed_genes.rds'))   # 5947 genes
gene.ensembl.id <- lapply(strsplit(gene.id, '[.]'), function(x){x[1]}) %>% unlist
# The gene ID is from GENCODE, which is not excatly the same as ENSEMBL ID. The part before [.] is the same as ENSEMBL ID.

library(biomaRt)
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl.37 <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = 37)

temp <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position','end_position'), 
              mart = ensembl, useCache = FALSE)
gene.mart <- temp[match(gene.ensembl.id, temp$ensembl_gene_id), ]  # Match ensembl id with chromosome positions.

# some gene id is from GRCh37, which has been deleted from CRCh38. We don't want to miss them.
left_gene <- gene.ensembl.id[is.na(gene.mart$ensembl_gene_id)]
temp.37 <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position','end_position'), 
                filters = 'ensembl_gene_id', values = left_gene, mart = ensembl.37, useCache = FALSE)
gene.mart[is.na(gene.mart$ensembl_gene_id), ] <- temp.37
gene.mart$chr <- as.factor(paste0('chr', gene.mart$chromosome_name))
saveRDS(gene.mart, file = paste0(processed_dir, "/gene_mart.rds"))

###########################
# 2. guide RNA positions
###########################
gRNA_indic_mat <- read.fst(paste0(processed_dir, "/gRNA_indicator_matrix.fst"))
gRNA.id <- colnames(gRNA_indic_mat)
temp.gRNA = strsplit(gRNA.id, ':')
gRNA.mart = data.frame(chr = lapply(temp.gRNA, function(x){x[1]}) %>% unlist)

temp.gRNA.p = do.call(rbind, lapply(temp.gRNA, function(x){unlist(strsplit(x[2], '-'))}))
gRNA.mart$start_position = as.numeric(temp.gRNA.p[, 1])
gRNA.mart$end_position = as.numeric(temp.gRNA.p[, 2])
gRNA.mart$gRNA_id = gRNA.id
gRNA.mart$mid_position = (gRNA.mart$start_position + gRNA.mart$end_position)/2  # Mid point of gRNA.
saveRDS(gRNA.mart, file = paste0(processed_dir, "/gRNA_mart.rds"))


###################################
# 3. gene-potential enhancer pairs
###################################

# From Gasperini paper, the distance between gRNA and gene is calculated as (TSS of gene - midpoint of gRNA). 
# They specifically select the gRNAs that are enhancers, that is, the position of gRNAs is ahead of the position of genes.
# The selection is performed chr by chr
chr.select = intersect(levels(gRNA.mart$chr), levels(gene.mart$chr))
# "chr1"  "chr10" "chr11" "chr12" "chr16" "chr17" "chr18" "chr19" "chr2"  
# "chr20" "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chrX" 
# 17 chromosomes

select.gRNA.gene.pair = NULL

for(chr in chr.select){
  gene.chr.id = which(gene.mart$chr == chr)
  gRNA.chr.id = which(gRNA.mart$chr == chr)
  gene.pos = gene.mart$start_position[gene.chr.id]
  gRNA.pos = gRNA.mart$mid_position[gRNA.chr.id]
  distance.temp = sapply(gene.pos, function(x){x - gRNA.pos})
  temp.id = which(distance.temp < 1000000 & distance.temp > 0, arr.ind = T)  
  # 1. Distance less than 1MB; 2. gRNA is in the front of gene region
  
  select.gRNA.gene.pair = rbind(select.gRNA.gene.pair, 
                                data.frame(gene.id = gene.id[gene.chr.id[temp.id[, 2]]], gRNA.id = gRNA.id[gRNA.chr.id[temp.id[, 1]]]))
  # cat(chr, ' is done! \n')
}
# dim(select.gRNA.gene.pair) # 1664 pairs of gene and gRNA
saveRDS(select.gRNA.gene.pair, file = paste0(processed_dir, "/select_gRNA_gene_pair.rds"))


####### Select transcription factor(TF) and their potiential enhancer 
##############################
# 4. select genes that are TF
##############################
tf.gene = read.csv(paste0(raw_data_dir, "/TF_human.csv"))
# preprocess TF_human.csv file. Delete lines that have the same ENSEMBL ID
tf.ensembl.id = unique(tf.gene$Animal.TFDB)
delete.tf.line = NULL
for(i in 1:length(tf.ensembl.id)){
  n = sum(tf.gene$Animal.TFDB == tf.ensembl.id[i])
  if(n > 1){
    cat(i, '\t')
    if(tf.ensembl.id[i] != ''){
      delete.tf.line = c(delete.tf.line, which(tf.gene$Animal.TFDB == tf.ensembl.id[i])[1])
    }
  }
}
# 24 	35 	183 	923 	1661 	
# delete tf.gene lines: 24, 193, 989, 1771

tf.gene.new = tf.gene[-delete.tf.line, ]
# We first select genes that are TF
match.temp = cbind(match(gene.mart$hgnc_symbol, tf.gene.new$Gene), match(gene.mart$ensembl_gene_id, tf.gene.new$Animal.TFDB))
# We match genes by their HGNC name and ENSEMBL name
match.temp = match.temp[rowSums(is.na(match.temp)) < 2, ]  # Remove genes that are not matched with HGNC name or ENSEMBL name

tf.gene.select = as.character(tf.gene$Gene[match.temp[, 1]])    # 406 genes
saveRDS(tf.gene.select, file = paste0(processed_dir, "/tf_gene_select.rds"))


######################################################
# 5. Find gRNAs that are pontential enhancers for TF 
######################################################
#temp <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position','end_position'), 
# mart = ensembl, useCache = FALSE)
### Find gRNAs that might be potential enhancers
tf.match.temp <- cbind(match(tf.gene.new$Gene, temp$hgnc_symbol), match(tf.gene.new$Animal.TFDB, temp$ensembl_gene_id))
tf.match <- tf.match.temp[, 2]
tf.match[is.na(tf.match.temp[, 2])] <- tf.match.temp[is.na(tf.match.temp[, 2]), 1]
undefine.id <- which(is.na(tf.match))
# [1]    6    8   17  202  203  204  205  214 1579 1631 1653 1666 1773
# Some of them are pseudogenes and do not make sense. We leave them first

tf.gene.mart <- temp[tf.match, ]
tf.gene.mart$chr <- paste0('chr', tf.gene.mart$chromosome_name)
saveRDS(tf.gene.mart, file = paste0(processed_dir, "/tf_gene_mart.rds"))

chr.select <- intersect(tf.gene.mart$chr, gRNA.mart$chr)  # 17 chromosome intersect
# [1] "chr16" "chr7"  "chr17" "chr19" "chr12" "chr2"  "chr8"  "chr20" "chr18" 
# "chr4"  "chrX"  "chr5"  "chr1"  "chr11" "chr10" "chr3"  "chr6" 
select.gRNA.tf.pair <- NULL

for(chr in chr.select){
  gene.chr.id <- which(tf.gene.mart$chr == chr)
  gRNA.chr.id <- which(gRNA.mart$chr == chr)
  gene.pos <- tf.gene.mart$start_position[gene.chr.id]
  gRNA.pos <- gRNA.mart$mid_position[gRNA.chr.id]
  distance.temp <- sapply(gene.pos, function(x){x - gRNA.pos})
  temp.id <- which(distance.temp < 1000000 & distance.temp > 0, arr.ind = T)  
  # 1. Distance less than 1MB; 2. gRNA is in the front of gene region
  
  select.gRNA.tf.pair <- rbind(select.gRNA.tf.pair, 
                              data.frame(ensembl.gene.id = tf.gene.mart$ensembl_gene_id[gene.chr.id[temp.id[, 2]]], 
                                         hgnc_symbol = tf.gene.mart$hgnc_symbol[gene.chr.id[temp.id[, 2]]], 
                                         gRNA.id = gRNA.id[gRNA.chr.id[temp.id[, 1]]]))
  cat(chr, ' is done! \n')
}
#dim(select.gRNA.tf.pair)
#[1] 366   3
# length(unique(select.gRNA.tf.pair$gRNA.id))
# includes 226 different gRNA
saveRDS(select.gRNA.tf.pair, file = paste0(processed_dir, "/select_gRNA_tf_pair.rds"))


##############################
# 6. negative control pairs
##############################
### gRNA that are not close to a transcription factor, pair with the set of genes on different chormosomes.

chr.select <- intersect(tf.gene.mart$chr, gRNA.mart$chr)  
num.neg.pair = NULL
neg.control.pair <- NULL
for(chr in chr.select){
  tf.chr.id <- which(tf.gene.mart$chr == chr)
  gRNA.chr.id <- which(gRNA.mart$chr == chr)
  tf.pos <- tf.gene.mart$start_position[tf.chr.id]
  gRNA.pos <- gRNA.mart$mid_position[gRNA.chr.id]
  distance.temp <- sapply(tf.pos, function(x){x - gRNA.pos})
  temp.gRNA = gRNA.id[gRNA.chr.id[which(rowSums(abs(distance.temp) > 1000000) == length(tf.pos))]]
  # Distance greater than 1MB
  temp.gene.ensembl = gene.ensembl.id[gene.mart$chr != chr]
  temp.hgnc = gene.mart$hgnc_symbol[gene.mart$chr != chr]
  
  neg.control.pair <- rbind(neg.control.pair, 
                            data.frame(gRNA.id = rep(temp.gRNA, length(temp.gene.ensembl)), 
                                       gene.ensembl.id = rep(temp.gene.ensembl, each = length(temp.gRNA)), 
                                       gene.hgnc.id = rep(temp.hgnc, each = length(temp.gRNA))))
  num.neg.pair = rbind(num.neg.pair, data.frame(chr = chr, gRNA = length(temp.gRNA), gene = length(temp.gene.ensembl)))
  #cat(chr, 'have',  length(temp.gRNA), 'gRNAs.', length(temp.gene.ensembl), ' genes not in ', chr, '. Done! \n')
}
saveRDS(neg.control.pair, file = paste0(processed_dir, '/neg_control_pair.rds'))
saveRDS(num.neg.pair, file = paste0(processed_dir, '/num_neg_pair.rds'))
