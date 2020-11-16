# Bulk RNA-seq
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
library(edgeR)

gene_names_in_use <- paste0(results_dir, "/all_results_with_names.fst") %>% read.fst() %>% filter(enh_names == "ARL15-enh") %>% pull(gene_names) %>% unique() %>% as.character()
bulk_data <- paste0(processed_dir, "/bulk_RNAseq.rds") %>% readRDS() %>% pluck("data", "arl15_enh") %>% select(-Chr, -Start, -End, -Strand, -Length) %>% filter(Geneid %in% gene_names_in_use)
bulk_info <- paste0(processed_dir, "/bulk_RNAseq.rds") %>% readRDS() %>% pluck("info", "arl15_enh")

group <- factor(bulk_info$region)
y <- DGEList(counts = select(bulk_data, -Geneid) %>% as.matrix, group = group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
genes_ordered <- topTags(qlf, n = nrow(bulk_data))

p_vals <- data.frame(gene_id = bulk_data$Geneid[row.names(genes_ordered) %>% as.integer()], p_value = genes_ordered$table$PValue, p_value_adj = genes_ordered$table$FDR)
saveRDS(object = p_vals, file = paste0(offsite_dir, "/results/xie/bulk_rna_seq/pvals_arl15_enh.rds"))
