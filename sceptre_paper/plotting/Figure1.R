args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper" else args[1]
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))

# QQ plot of Gasperini NTC p-values
ci = 0.95
subsampling_factor = 250
p = original_results %>%
  rename(pvalue = pvalue.raw) %>%
  filter(site_type == "NTC", beta < 0) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  # mutate(pvalue = ifelse(pvalue < 1e-50, 0, pvalue)) %>%
  filter(-log10(expected) > 2 | row_number() %% subsampling_factor == 0) %>%
  ggplot(aes(x = expected, y = pvalue, ymin = clower, ymax = cupper)) + 
  geom_point(aes(colour = site_type), size = 3, alpha = 0.5) + 
  geom_ribbon(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("dodgerblue")) + 
  xlab("Expected null p-value") + 
  ylab("Observed p-value for gene-gRNA pair") + 
  scale_x_continuous(trans = revlog_trans(base = 10)) + 
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  theme_bw() + theme(text = element_text(size = 14),
                     legend.position = "none",
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)

# QQ plot of Xie ARL15 p-values


# Dispersion estimation

KS_pvals = original_results %>% 
  rename(gene_id = ENSG) %>%
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
p = df %>% arrange(desc(KS_pval)) %>% mutate(KS_pval = ifelse(KS_pval < 1e-5, 1e-5, KS_pval)) %>% 
  ggplot() + geom_point(aes(x = mu, y = disp_raw, colour = KS_pval)) + 
  geom_point(aes(x = mu, y = disp_raw), size = 4, shape = 21, data = df %>% filter(KS_pval < 0.05/num_genes)) + 
  geom_line(aes(x = mu, y = disp_coeffs[1] + disp_coeffs[2]/mu), colour = "black", linetype = "dashed") + 
  scale_colour_continuous(trans = "log10", name = "p-value miscalibration,\nper gene") +
  scale_x_log10() + scale_y_log10() + 
  xlab("Mean gene expression") + ylab("Dispersion of gene expression") + 
  guides(colour = guide_colorbar(title.position = "top", direction = "horizontal", order = 1)) +
  guides(shape = guide_legend(title = NULL, direction = "vertical", order = 2)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = c(0.65, 0.175),
                     panel.border = element_blank(), axis.line = element_line(),
                     legend.box = "horizontal",  legend.key = element_rect(colour = "transparent", fill = "white"),
                     text = element_text(size = 14),
                     legend.title = element_text(size = 12),
                     legend.background = element_rect(fill = "transparent", colour = NA))
plot(p)

# confounding for Xie
covariates_xie$guide_count = rowSums(grna_indicator_matrix_xie)
p = covariates_xie %>% 
  mutate(total_umis = 10^(log_n_umis)) %>%
  ggplot(aes(x = total_umis, y = guide_count)) +
  stat_density2d(aes(fill = ..density..^0.25), adjust = 0.5, geom = "raster", contour = FALSE) +
  # geom_point(data = covariates_xie %>% filter(computed_density < 0.001)) +
  # stat_density_2d(aes(fill = stat(log(level))), geom = 'polygon', bins = 20) +
  # scale_fill_viridis_c(name = "density") +
  geom_smooth(method = "lm", colour = "black", se = FALSE, linetype = "dashed") + 
  scale_x_log10(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  xlab("Total UMIs per cell") + ylab("Total gRNAs per cell") + 
  scale_fill_gradient(low = "white", high = "dodgerblue4") +
  theme_bw() + theme(legend.position = "none", text = element_text(size = 14), 
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)

# confounding for Gasperini
