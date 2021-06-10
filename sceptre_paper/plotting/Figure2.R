# Reproduce Figure 2 from Katsevich, Barry, and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
offsite_dir <- if (is.na(args[2])) "/Volumes/tims_new_drive/research/sceptre_files" else args[2]
fig2_dir <- paste0(code_dir, "/sceptre_paper/manuscript/figures/Figure2")
require(sceptre)
require(katsevich2020)

# Select the gene and gRNA to use in the plot
gene_to_plot <- "ENSG00000135390"
grna_to_plot <- "chr12.1893_top_two"
sceptre_function_args <- paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_gasp/sceptre_function_args.R")
p <- plot_skew_t_gene_gRNA(gene_id = gene_to_plot, gRNA_id = grna_to_plot, sceptre_function_args = sceptre_function_args, offsite_dir = offsite_dir)
plot(p$plot)
ggsave(filename = paste0(fig2_dir, "/Figure2c.png"), plot = p$plot, device = "png", width = 3.5, height = 2.5, dpi = "retina")
