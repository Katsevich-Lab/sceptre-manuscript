args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
suppressPackageStartupMessages(library(Gmedian))
suppressPackageStartupMessages(library(gridExtra))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_gasp/sceptre_function_args.R" else args[2]
raw_data_dir <- paste0(offsite_dir, "/data/gasperini/raw")
sim_data_dir <- paste0(offsite_dir, "/data/simulations")
source(param_file)

# We select a "representative" Gaspeini gene on which to perform the analysis. This gene is the geometric median of the (log(geometric_mean), log(size)) scatter plot.
file_names <- list.files(gene_precomp_dir)
gene_size_unreg_files <- paste0(gene_precomp_dir, "/", grep(pattern = 'gene_size_unreg_[0-9]+.rds', x = file_names, value = TRUE))
gene_geom_mean_files <- paste0(gene_precomp_dir, "/", grep(pattern = 'gene_geom_mean_[0-9]+.rds', x = file_names, value = TRUE))
load_distributed_vector <- function(f_names) f_names %>% map(readRDS) %>% reduce(c)
theta <- load_distributed_vector(gene_size_unreg_files)
genes_log_gmean <- load_distributed_vector(gene_geom_mean_files)
all(colnames(theta) == colnames(genes_log_gmean)) # make sure genes are in same order in both vectors
mean_dispersion_df <- tibble(gene_id = names(genes_log_gmean), log_geom_mean = set_names(genes_log_gmean, NULL), log_size = log(set_names(theta, NULL)))
g_median <- Gmedian(select(mean_dispersion_df, -gene_id))
row_no <- which.min(Mod((g_median[1] + g_median[2] * 1i) - (mean_dispersion_df$log_geom_mean + mean_dispersion_df$log_size * 1i)))
my_gene <- slice(mean_dispersion_df, row_no) %>% pull(gene_id)
p1 <- ggplot(data = mean_dispersion_df, mapping = aes(x = log_geom_mean, y = log_size)) + geom_point(alpha = 0.3) + theme_bw() + xlab("log(geometric mean)") + ylab("log(size)") + geom_point(inherit.aes = FALSE, data = slice(mean_dispersion_df, row_no), mapping = aes(x = log_geom_mean, y = log_size), color = "red", size = 2) + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + annotate(geom = "text", label = my_gene, col = "red", x = -0.2, y = 4.5) + ggtitle("All genes")

# We similarly select a "representative" NTC gRNA on which to perform the analysis. This will be the NTC gRNA with median expression.
gRNA_site_types <- suppressWarnings(read_tsv(paste0(raw_data_dir, "/GSE120861_all_deg_results.at_scale.txt"), col_types = "cddddddccccciiciiccc")) %>% select(gRNA_group, site_type) %>% distinct()
ntc_ids <- filter(gRNA_site_types, site_type == "NTC") %>% pull(gRNA_group)
rm(gRNA_site_types)
ntc_indicator_mat <- read.fst(path = gRNA_indicator_matrix_fp, columns = ntc_ids)
ntc_means <- apply(ntc_indicator_mat, 2, mean)
my_gRNA <- names(which.min(abs(ntc_means - median(ntc_means))))
p2 <- ggplot(data = tibble(ntc_means), mapping = aes(x = 100 * ntc_means)) + geom_histogram(bins = 8, color="black", fill="grey") + theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("Percent cells infected") + geom_vline(xintercept = 100 * ntc_means[[my_gRNA]], col = "red") + ylab("") + annotate(geom = "text", x = .65, y = 11, label = my_gRNA, col = "red") + ggtitle("NTC gRNAs")

p_save <- grid.arrange(p1, p2, nrow = 1)
ggsave(filename = paste0(offsite_dir, "/figures/simulation_chosen_gRNA_gene.pdf"), plot = p_save, device = "pdf", scale = 1)

to_save <- c(gene = my_gene, gRNA = my_gRNA)
saveRDS(object = to_save, file = paste0(sim_data_dir, "/simulation_gRNA_gene.rds"))
