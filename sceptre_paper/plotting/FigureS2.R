# Reproduce Figure S2 from Katsevich, Barry, and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
require(katsevich2020)
require(cowplot)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
figS2_dir <- paste0(manuscript_figure_dir, "/FigureS2")

# Collate candidate enhancer p-values from negative binomial and
# conditional randomization approaches
df = likelihood_results_gasp %>%
  filter(site_type == "DHS") %>%
  select(gene_id, grna_group, pvalue) %>%
  dplyr::rename(old_pvalue = pvalue) %>%
  left_join(resampling_results_gasp %>%
              filter(site_type == "DHS") %>%
              select(gene_id, grna_group, p_value, xi, alpha, omega, nu) %>%
              dplyr::rename(new_pvalue = p_value),
            by = c("grna_group", "gene_id"))

# First, locate the old and new p-value that are closest to one another. Consider a subset of pairs chosen randomly.
df_samp <- df %>% filter(!(gene_id %in% c("ENSG00000109654", "ENSG00000099860", "ENSG00000223609")), grna_group != "chr6.1172_top_two") %>% mutate(p_dif = new_pvalue - old_pvalue) # sample_n(tbl = df, size = 100000, replace = FALSE) %>% mutate(p_dif = new_pvalue - old_pvalue)

same_p_val <- filter(df_samp, new_pvalue >= 0.05 & new_pvalue <= 0.3) %>% mutate(d_from_z_squared = (xi)^2 + (alpha)^2 + (omega-1)^2 + (1/nu)^2) %>% summarize(min_idx = which.min(d_from_z_squared), gene_id = gene_id[min_idx], gRNA_id = grna_group[min_idx], old_pvalue = old_pvalue[min_idx], new_pvalue = new_pvalue[min_idx])
sceptre_bigger <- filter(df_samp, new_pvalue >= 0.05 & new_pvalue <= 0.3) %>% summarize(min_idx = which.max(p_dif), gene_id = gene_id[min_idx], gRNA_id = grna_group[min_idx], old_pvalue = old_pvalue[min_idx], new_pvalue = new_pvalue[min_idx])
sceptre_smaller <- filter(df_samp, old_pvalue >= 0.05 & old_pvalue <= 0.4) %>% summarize(min_idx = which.min(p_dif), gene_id = gene_id[min_idx], gRNA_id = grna_group[min_idx], old_pvalue = old_pvalue[min_idx], new_pvalue = new_pvalue[min_idx])

# three gene/gRNA pairs to illustrate the differences in the two approaches
genes_to_plot <- c(same_p_val$gene_id, sceptre_smaller$gene_id, sceptre_bigger$gene_id)
grnas_to_plot <- c(same_p_val$gRNA_id, sceptre_smaller$gRNA_id, sceptre_bigger$gRNA_id)

# panel a: comparing the two ways of computing p-values
p0 = df %>%
  mutate(highlighted =
           (gene_id == genes_to_plot[1] & grna_group == grnas_to_plot[1]) |
           (gene_id == genes_to_plot[2] & grna_group == grnas_to_plot[2]) |
           (gene_id == genes_to_plot[3] & grna_group == grnas_to_plot[3])) %>%
  arrange(highlighted) %>%
  filter(old_pvalue > 1e-10, new_pvalue > 1e-10, !is.na(alpha), alpha < 4) %>%
  ggplot(aes(x = old_pvalue, y = new_pvalue, colour = highlighted)) +
  geom_point(alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  xlab("Hafemeister NB p-value") + ylab("SCEPTRE p-value") +
  scale_colour_manual(values = c("darkslategray4", "firebrick3")) +
  theme_bw() + theme(legend.position = "none",
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line())
plot(p0)

# panel b: NB p-value more significant
p1 <- plot_skew_t_gene_gRNA(genes_to_plot[3], grnas_to_plot[3], interval = c(-6,6))
p1_plot <- p1$plot + ggtitle("\n") + theme(legend.position = "none")

# panel c: SCEPTRE p-value more significant
p2 <- plot_skew_t_gene_gRNA(genes_to_plot[2], grnas_to_plot[2], c(-4, 4))
p2_plot <- p2$plot + xlim(c(-4,4)) + ggtitle("\n") + theme(legend.position = c(.85, .8), legend.background = element_rect(color = "transparent", fill = "transparent"))

# panel d: two p-values about equal
p3 <- plot_skew_t_gene_gRNA(genes_to_plot[1], grnas_to_plot[1], c(-4,4))
p3_plot <- p3$plot + theme(legend.position = "none") + ggtitle("\n")

# combine the panels
left_col <- plot_grid(p0, p1_plot, align = "h", ncol = 1, labels = c("a", "c"))
right_col <- plot_grid (p2_plot, p3_plot, align = "vh", ncol = 1, labels = c("b", "d"))
plot_final <- plot_grid(left_col, right_col, ncol = 2)
ggsave(filename = paste0(figS2_dir, "/FigureS2a-d.pdf"), plot = plot_final, device = "pdf", scale = 1, width = 8.0, height = 7)
