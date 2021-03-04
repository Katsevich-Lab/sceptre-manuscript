# Reproduce Figure 3 from Katsevich, Barry, and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
require(katsevich2020)
require(cowplot)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))

resampling_results_xie <- resampling_results_xie %>% filter(enh_names == "MYB-enh-3")
resampling_results_xie_for_bulk <- resampling_results_xie %>% select(p_value, gene_names) %>% mutate(p_value_adj = p.adjust(p_value, method = "BH"), rejected_sceptre = p_value_adj < 0.1)
p_vals_bulk <- mutate(p_vals_bulk_myb3_enh3, rejected_bulk = p_value_adj < 0.1)

# resampling_results_xie <- resampling_results_xie %>% filter(enh_names == "ARL15-enh")
# resampling_results_xie_for_bulk <- resampling_results_xie %>% select(p_value, gene_names) %>% mutate(p_value_adj = p.adjust(p_value, method = "BH"), rejected_sceptre = p_value_adj < 0.1)
# p_vals_bulk <- mutate(p_vals_bulk, rejected_bulk = p_value_adj < 0.1)

to_plot <- inner_join(p_vals_bulk, resampling_results_xie_for_bulk, by = "gene_names") %>% rename(bulk_pval_adj = p_value_adj.x, bulk_pval = p_value.x, sceptre_pval_adj = p_value_adj.y, sceptre_pval = p_value.y)

table(to_plot$rejected_sceptre, to_plot$rejected_bulk)
p_e <- ggplot(data = to_plot, mapping = aes(x = bulk_pval, y = sceptre_pval)) +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = c("grey60", "firebrick3")) +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  geom_vline(xintercept = filter(to_plot, !rejected_bulk) %>% pull(bulk_pval) %>% min(), col = "royalblue4", linetype = "dashed") +
  geom_hline(yintercept = filter(to_plot, !rejected_sceptre) %>% pull(sceptre_pval) %>% min(), col = "royalblue4", linetype = "dashed") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.background = element_rect(fill = "transparent", color = NA),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(), legend.position = "none") +
  xlab("Bulk RNA-seq p-value") + ylab("SCEPTRE p-value") + ggtitle("Bulk RNA-seq validation for MYB-enh3")
