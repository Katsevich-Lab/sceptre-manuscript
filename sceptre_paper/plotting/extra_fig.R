args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/sceptre-manuscript/SCEPTRE/" else args[1]
require(katsevich2020)
require(cowplot)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))

resampling_results_xie_myb <- resampling_results_xie %>% filter(enh_names == "MYB-enh-3")
# add colum for two-sided p-value
get_two_sided_p <- Vectorize(function(skew_t_success, z_value, xi, omega, alpha, nu) {
  if (skew_t_success) {
    dp <- c(xi, omega, alpha, nu)
    pmax(.Machine$double.eps, pst(x = -abs(z_value), dp = dp, method = 2, rel.tol = .Machine$double.eps) +
          (1 - pst(x = abs(z_value), dp = dp, method = 2, rel.tol = .Machine$double.eps)))
  } else NA
})

resampling_results_xie_myb_fixed <- resampling_results_xie_myb %>% mutate(two_tailed_p = get_two_sided_p(skew_t_fit_success, z_value, xi, omega, alpha, nu)) %>%
  select(p_value = two_tailed_p, gene_names) %>% mutate(p_value_adj = p.adjust(p_value, method = "BH"), rejected_sceptre = p_value_adj < 0.1)
p_vals_bulk <- mutate(p_vals_bulk_myb3_enh3, rejected_bulk = p_value_adj < 0.1)
to_plot <- inner_join(p_vals_bulk, resampling_results_xie_myb_fixed, by = "gene_names") %>% rename(bulk_pval_adj = p_value_adj.x, bulk_pval = p_value.x, sceptre_pval_adj = p_value_adj.y, sceptre_pval = p_value.y)

table(to_plot$rejected_sceptre, to_plot$rejected_bulk)
p_a <- ggplot(data = to_plot, mapping = aes(x = bulk_pval, y = sceptre_pval)) +
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


# qq-plot
df1 <- to_plot %>% pivot_longer(cols = c("bulk_pval", "sceptre_pval"), names_to = "method", values_to = "pvalue") %>%
  group_by(method) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r],
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>% mutate()
p_b <- df1 %>%
  ggplot(aes(x = expected, y = pvalue, col = method)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_x_continuous(trans = revlog_trans(base = 10)) + scale_y_continuous(trans = revlog_trans(base = 10)) +
  xlab(expression(paste("Expected null p-value"))) +
  ylab(expression(paste("Observed p-value"))) +
  ggtitle("SCEPTRE vs. bulk p-values for MYB-enh3") +
  theme_bw() + theme(legend.background = element_rect(fill = "transparent", colour = NA),
                     legend.title = element_blank(),
                     panel.grid = element_blank(),
                     strip.background = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
