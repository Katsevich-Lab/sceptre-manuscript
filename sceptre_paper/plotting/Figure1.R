# Reproduce Figure 1 from Katsevich, Barry, and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
require(katsevich2020)
require(cowplot)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
fig1_dir <- paste0(manuscript_figure_dir, "/Figure1")

# Subfigure a: Gasperini NTC p-values and Xie negative control p-value qqplots
p_thresh <- 1e-08
temp_down = rbind(  # downsample Gasperini NTC p-values so that Xie and Gasperini have the same number of NTC pairs.
  original_results_gasp %>%
    rename(pvalue = pvalue.raw) %>%
    filter(site_type == "NTC", beta < 0) %>% sample_n(85000) %>%
    mutate(r = rank(pvalue), expected = ppoints(n())[r],
           clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
           cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
    filter(-log10(expected) > 0 | row_number() %% subsampling_factor == subsampling_factor) %>%
    mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>%
    select(pvalue, expected, clower, cupper) %>%
    mutate(data = 'Gasperini'),
  original_results_xie %>%
    rename(pvalue = raw_p_val) %>%
    filter(site_type == "negative_control") %>%
    mutate(r = rank(pvalue), expected = ppoints(n())[r],
           clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
           cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
    filter(-log10(expected) > 0 | row_number() %% subsampling_factor == subsampling_factor) %>%
    mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>% 
    select(pvalue, expected, clower, cupper) %>%
    mutate(data = 'Xie'))

p_a = ggplot(data = temp_down, aes(x = expected, y = pvalue, ymin = clower, ymax = cupper, color = data)) +
  geom_point(size = 1.2, alpha = 0.5) + 
  scale_color_manual(values = c(plot_colors[["gasperini_nb"]], plot_colors[["hypergeometric"]])) + 
  geom_ribbon(alpha = 0.5, color = NA, fill = 'gray') +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Expected null p-value") +
  ylab("Observed p-value") +
  ggtitle("P-value miscalibration") +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  theme_bw() + theme(
    #legend.position = "none",
    legend.position = c(0.22, 0.77),
    legend.title = element_blank(), 
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.background = element_blank(),
    legend.margin = margin(-0.25,0,0,0, unit="cm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

hjust <- -0.5
vjust <- 3
# Subfigure b
gasp_confounding_df <- covariates_gasp %>% summarize(total_umis = ifelse(is.na(lg_total_umis), 0, exp(lg_total_umis)), guide_count = ifelse(is.na(lg_guide_count), 0, exp(lg_guide_count)))
rho <- cor(log(gasp_confounding_df$total_umis), gasp_confounding_df$guide_count) %>% signif(2)
cor.test(log(gasp_confounding_df$total_umis), gasp_confounding_df$guide_count)

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
cor.test(xie_confounding_df$total_umis, xie_confounding_df$guide_count)
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
p_top <- plot_grid(p_a, p_b, p_c, labels = c("a", "b", "c"), nrow = 1)
p_bottom <- plot_grid(p_d, p_e, rel_widths = c(1,2), nrow = 1, labels = c("d", "e"))
p_final <- plot_grid(p_top, p_bottom, nrow = 2, rel_heights = c(0.8, 1), align = "h")
ggsave(filename = paste0(fig1_dir, "/figure1_abce.pdf"), plot = p_final, device = "pdf", scale = 1, width = 9.5, height = 7)
