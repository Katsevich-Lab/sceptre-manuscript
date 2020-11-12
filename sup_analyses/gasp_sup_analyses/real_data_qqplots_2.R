# Load the gasperini, (HF) negative binomial, and sceptre results

args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
require(katsevich2020)
source(paste0(code_dir, "analysis_drivers_gasp/file_paths_to_dirs.R"))
source(paste0(code_dir, "sup_analyses/xie_sup_analyses/aux_objects.R"))

gasperini_results <- paste0(processed_dir, "/original_results.fst") %>% read.fst()
sceptre_results <- paste0(results_dir, "/resampling_results.fst") %>% read.fst()
negbin_results <- paste0(results_dir_negative_binomial, "/likelihood_results.fst") %>% read.fst()
gene_orig_results <- paste0(offsite_dir, "/results/comparison_to_gene/resampling_results.tsv") %>% read_tsv()

df_NTC <- rbind(select(gasperini_results, gene_id, grna_group, pvalue = pvalue.raw, site_type) %>% mutate(method = "Original NB"),
                select(sceptre_results, gene_id, grna_group, pvalue = p_value, site_type) %>% mutate(method = "SCEPTRE"),
                select(negbin_results,  gene_id, grna_group, pvalue, site_type) %>% mutate(method = "Regularized NB")) %>% filter(site_type == "NTC") %>% mutate_at(.vars = c("gene_id", "grna_group", "site_type", "method"), .funs = factor)

# compute information for overall QQ plot
ci <- 0.95
subsampling_factor <- 250

df1 <- df_NTC %>%
  group_by(method) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r],
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup()

# compute information for QQ plot by gRNA
df2 <- df_NTC %>%
  group_by(grna_group, method) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r],
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>%
  group_by(r, expected, clower, cupper, method) %>%
  summarise(pmedian = median(pvalue), plower = quantile(pvalue, 0.025), pupper = quantile(pvalue, 0.975)) %>%
  ungroup()

# Overall calibration on real data
p1 <- df1 %>% filter(-log10(expected) > 2 | row_number() %% subsampling_factor == 0) %>% mutate(clower = ifelse(method == "SCEPTRE", clower, NA),
                                                       cupper = ifelse(method == "SCEPTRE", cupper, NA)) %>%
  mutate(pvalue = ifelse(pvalue < 1e-8, 1e-8, pvalue)) %>%
  ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) +
  geom_point(aes(color = method), size = 1, alpha = 0.5) +
  geom_ribbon(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c(plot_colors[["original_nb"]], plot_colors[["improved_nb"]], plot_colors[["sceptre"]])) +
  scale_x_continuous(trans = revlog_trans(base = 10)) + scale_y_continuous(trans = revlog_trans(base = 10)) +
  xlab(expression(paste("Expected null p-value"))) +
  ylab(expression(paste("Observed p-value"))) +
  ggtitle("All real gene-NTC pairs") +
  theme_bw() + theme(legend.position = c(0.75,0.20),
                     legend.background = element_rect(fill = "transparent", colour = NA),
                     legend.text = element_text(size = 11),
                     legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5, size = 12),
                     panel.grid = element_blank(),
                     strip.background = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave(filename = paste0(offsite_dir, "/figures/gasp_qqplot_combined.pdf"), plot = p1, device = "pdf", scale = 1, width = 5, height = 4)

# Calibration per gene on real data
p2 <- df2 %>%
  mutate_at(c("expected", "pmedian", "plower", "pupper"), ~ifelse(.<1e-6, 1e-6,.)) %>%
  ggplot(aes(x = expected, y = pmedian)) +
  geom_line(aes(colour = method), size = 1) +
  geom_ribbon(aes(ymin = plower, ymax = pupper, fill = method), alpha = 0.4) +
  geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c(plot_colors[["original_nb"]], plot_colors[["improved_nb"]], plot_colors[["sceptre"]])) +
  scale_fill_manual(values = c(plot_colors[["original_nb"]], plot_colors[["improved_nb"]], plot_colors[["sceptre"]])) +
  facet_wrap(method~., scales = "fixed", nrow = 1) +
  scale_x_continuous(trans = revlog_trans(base = 10), breaks = c(1,1e-2,1e-4)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  ggtitle("Real gene-NTC pairs for each NTC") +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(hjust = 0.5, size = 12),
                     panel.spacing.x = unit(1.25, "lines"),
                     axis.title = element_blank(),
                     strip.background = element_blank(),
                     strip.text = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line())
ggsave(filename = paste0(offsite_dir, "/figures/gasp_qqplot_per_gene.pdf"), plot = p2, device = "pdf", scale = 1, width = 6, height = 3)

#################################
joined <- left_join(x = gene_orig_results %>% filter(method == "conditional_randomization") %>% select(pval_old = corrected_pvalue_st, pair_id, quality_rank_grna, site_type), y = sceptre_results %>% select(pval_new = p_value, pair_id), by = "pair_id")

joined_adj <- joined %>% filter(quality_rank_grna == "top_two" & site_type == "DHS") %>% mutate(pval_old_adj = p.adjust(pval_old, method = "BH"), pval_new_adj = p.adjust(pval_new, method = "BH"), reject_old = pval_old_adj <= 0.1, reject_new = pval_new_adj <= 0.1) 

p_thresh_old <- filter(joined_adj, reject_old) %>% summarize(max(pval_old)) %>% pull()
p_thresh_new <- filter(joined_adj, reject_new) %>% summarize(max(pval_new)) %>% pull()

p <- ggplot(data = joined %>% filter(pval_old <= 0.1 | pval_new <= 0.1 | row_number() %% subsampling_factor == 0), mapping = aes(x = pval_old, y = pval_new)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_hline(yintercept = p_thresh_new, col = "blue") + geom_vline(xintercept = p_thresh_old, col = "blue")
