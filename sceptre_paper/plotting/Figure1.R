# Reproduce Figure 1 from Katsevich, Barry, and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper" else args[1]
require(katsevich2020)
require(ggpubr)
require(cowplot)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
fig1_dir <- paste0(manuscript_figure_dir, "/Figure1")

# Subfigure a: Gasperini NTC p-values and Xie ARL15 p-value qqplots
p_thresh <- 1e-21
pa.1 <- original_results_gasp %>%
  rename(pvalue = pvalue.raw) %>%
  filter(site_type == "NTC", beta < 0) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r],
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  filter(-log10(expected) > 2 | row_number() %% subsampling_factor == subsampling_factor) %>%
  mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>%
  ggplot(aes(x = expected, y = pvalue, ymin = clower, ymax = cupper)) +
  geom_point(col = plot_colors[["gasperini_nb"]], size = 2, alpha = 0.5) +
  geom_ribbon(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("") +
  ylab("") +
  ggtitle("Gasperini") +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  theme_bw() + theme(
                     legend.position = "none",
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     plot.title = element_text(hjust = 0.5))

p_vals_hypergeo <- original_results_xie %>% pluck("arl15_enh")
arl15_id <- filter(resampling_results_xie, gene_names == "ARL15") %>% pull(gene_id) %>% as.character()
p_vals_hypergeo <- p_vals_hypergeo[names(p_vals_hypergeo) %in% resampling_results_xie$gene_id]
p_vals_hypergeo <- p_vals_hypergeo[!(names(p_vals_hypergeo) %in% arl15_id)]

pa.2 <- tibble(pvalue = p_vals_hypergeo, gene_id = names(p_vals_hypergeo)) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r],
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ggplot(aes(x = expected, y = pvalue, ymin = clower, ymax = cupper)) +
  geom_point(col = plot_colors[["hypergeometric"]], size = 1, alpha = 0.5) +
  geom_ribbon(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Expected null p-value") +
  ylab("") +
  ggtitle("Xie") +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  theme_bw() + theme(
                     legend.position = "none",
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     plot.title = element_text(hjust = 0.5))
top_left <- ggarrange(pa.1, pa.2, ncol = 1)
pa_annotated <- annotate_figure(top_left, left = text_grob("Observed p-value gene-gRNA pair", rot = 90, vjust = 2.5))

hjust <- -0.5
vjust <- 3
# Subfigure b
gasp_confounding_df <- covariates_gasp %>% summarize(total_umis = ifelse(is.na(lg_total_umis), 0, exp(lg_total_umis)), guide_count = ifelse(is.na(lg_guide_count), 0, exp(lg_guide_count)))
rho <- cor(log(gasp_confounding_df$total_umis), gasp_confounding_df$guide_count) %>% signif(2)
text_to_add <- paste0("~ rho == ", rho)

p_b <- gasp_confounding_df %>% filter(guide_count > 0, row_number() %% 10 == 0) %>%
  ggplot(aes(x = total_umis, y = guide_count)) +
  stat_density2d(aes(fill = ..density..^0.25), adjust = 0.5, geom = "raster", contour = FALSE) +
  geom_smooth(formula = y~x, method = "lm", colour = "black", se = FALSE, linetype = "dashed") +
  scale_x_log10(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Gasperini confounding") +
  xlab("Total UMIs per cell") + ylab("Total gRNAs per cell") +
  scale_fill_gradient(low = "white", high = plot_colors[["gasperini_nb"]]) +
  theme_bw() + theme(legend.position = "none",
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     plot.title = element_text(hjust = 0.5)) +
  annotate(geom = "text", x = 0, y = Inf, label = as.character(text_to_add), parse = TRUE, size = 5, hjust = hjust, vjust = vjust)

# subfigure c
guide_count <- rowSums(grna_indicator_matrix_xie)
guide_count[is.na(guide_count)] <- 0
umi_count <- exp(covariates_xie$log_n_umis)
xie_confounding_df <- tibble(total_umis = umi_count, guide_count = guide_count) %>% filter(guide_count <= 50)
rho <- cor(xie_confounding_df)[1,2] %>% signif(2)
text_to_add <- paste0("~ rho == ", rho)

p_c <- xie_confounding_df %>% filter(guide_count > 0, row_number() %% 10 == 0) %>%
  ggplot(aes(x = total_umis, y = guide_count)) +
  stat_density2d(aes(fill = ..density..^0.25), adjust = 0.5, geom = "raster", contour = FALSE) +
  geom_smooth(formula = y~x, method = "lm", colour = "black", se = FALSE, linetype = "dashed") +
  scale_x_log10(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Total UMIs per cell") + ylab("") +
  ggtitle("Xie confounding") +
  scale_fill_gradient(low = "white", high = plot_colors[["hypergeometric"]]) +
  theme_bw() + theme(legend.position = "none",
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     axis.line = element_line()) + annotate(geom = "text", x = 0, y = Inf, label = as.character(text_to_add), parse = TRUE, size = 5, hjust = hjust, vjust = vjust)

# subfigure d
p_d <- ggplot() + theme_bw() + ggtitle("Confounding structure") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

# subfigure e
KS_pvals = original_results_gasp %>%
  filter(site_type == "NTC") %>%
  group_by(gene_id) %>%
  summarise(KS_pval = ks.test(pvalue.raw, "punif")$p.value) %>%
  ungroup()
df = disp_table %>%
  rename(disp_raw = disp) %>%
  mutate(disp_shrunk = disp_coeffs[1] + disp_coeffs[2]/mu) %>%
  inner_join(KS_pvals, by = "gene_id") %>%
  mutate(KS_pval_adj = p.adjust(KS_pval, "fdr"))
num_genes = nrow(df)
p_e <- df %>% arrange(desc(KS_pval)) %>% mutate(KS_pval = ifelse(KS_pval < 1e-5, 1e-5, KS_pval)) %>%
  ggplot() + geom_point(aes(x = mu, y = disp_raw, colour = KS_pval)) +
  geom_point(aes(x = mu, y = disp_raw), size = 4, shape = 21, data = df %>% filter(KS_pval < 0.05/num_genes)) +
  geom_line(aes(x = mu, y = disp_coeffs[1] + disp_coeffs[2]/mu), colour = "black", linetype = "dashed") +
  scale_colour_continuous(trans = "log10", name = "p-value miscalibration, per gene") +
  scale_x_log10() + scale_y_log10() +
  xlab("Mean gene expression") + ylab("Dispersion of gene expression") + ggtitle("Dispersion estimation (Gasperini data)") +
  guides(colour = guide_colorbar(title.position = "top", direction = "horizontal", order = 1)) +
  guides(shape = guide_legend(title = NULL, direction = "vertical", order = 2)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = c(0.7, .15),
                     panel.border = element_blank(), axis.line = element_line(),
                     legend.box = "horizontal",  legend.key = element_rect(colour = "transparent", fill = "white"),
                     plot.title = element_text(hjust = 0.5),
                     legend.background = element_rect(fill = "transparent", colour = NA))

# Combine everything through cowplot
p_top <- plot_grid(pa_annotated, p_b, p_c, labels = c("a", "b", "c"), nrow = 1)
p_bottom <- plot_grid(p_d, p_e, rel_widths = c(1,2), nrow = 1, labels = c("d", "e"))
p_final <- plot_grid(p_top, p_bottom, nrow = 2, rel_heights = c(0.8, 1), align = "h")
ggsave(filename = paste0(fig1_dir, "/figure1_abce.pdf"), plot = p_final, device = "pdf", scale = 1, width = 9.5, height = 7)
