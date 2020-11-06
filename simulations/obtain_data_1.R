args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_xie/sceptre_function_args.R" else args[2]
source(param_file)

if (FALSE) {
# check that the number of cells in the cell-gene expression matrix matches the number of cells in the gRNA indicator matrix. Also, check that the gRNA indicator matrix is free of NAs when we subset to cell_subset.
gRNA_indic_mat <- read.fst(gRNA_indicator_matrix_fp)
nrow(gRNA_indic_mat) == nrow(cell_gene_expression_matrix)
gRNA_indic_mat[cell_subset,] %>% is.na() %>% sum() 

randomly_select_in_quantile <- function(vect) {
  percent_express_quants <- quantile(vect, c(0.7,0.95))
  possible_vals <- names(which(vect <= percent_express_quants[["95%"]] & vect >= percent_express_quants[["70%"]]))
  sample(possible_vals, 1)
}

# survey the genes; pick one with high expression (but not too high) for the example
set.seed(1)

gene_expression_p <- big_apply(cell_gene_expression_matrix, function(X, ind) colMeans(X[,ind] >= 1)) %>% unlist() * 100
names(gene_expression_p) <- ordered_gene_ids
my_gene <- randomly_select_in_quantile(gene_expression_p)
expressions <- cell_gene_expression_matrix[,which(my_gene == ordered_gene_ids)]

gRNA_percent_express <- apply(gRNA_indic_mat, 2, function(col) mean(col, na.rm = TRUE)) * 100
my_gRNA <- randomly_select_in_quantile(gRNA_percent_express)
gRNA_indics <- gRNA_indic_mat[,my_gRNA] %>% na.omit()
bulk_region_names <- readRDS(paste0(processed_dir, "/bulk_region_names.rds"))
my_gRNA %in% bulk_region_names
}

my_gRNA <- "chr18:9705278-9705678"
my_gene <- "ENSG00000109381.19"

# save the selected gene and gRNA
saveRDS(object = c(gRNA = my_gRNA, gene = my_gene), paste0(offsite_dir, "/data/simulations/simulation_gRNA_gene.rds"))

# Note: our gene is on chromosome 4 and our gRNA is on chromosome 18; moreover, our gRNA does not target either ARL15 or MYB.