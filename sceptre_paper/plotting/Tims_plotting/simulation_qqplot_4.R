args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
offsite_dir <- if (is.na(args[2])) "/Volumes/tims_new_drive/research/sceptre_files" else args[2]
require(katsevich2020)
source(paste0(code_dir, "/sceptre_paper/plotting/aux_objects.R"))

truncate_thresh <- 1e-7
qq_data <- simulation_results %>%
  group_by(method, dataset_id) %>%
  mutate(r = rank(p_value), expected = ppoints(n())[r],
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>%
  mutate(pvalue = ifelse(p_value < truncate_thresh, truncate_thresh, p_value),
         facet_label = factor(x = as.character(dataset_id), levels = c("2", "1", "3", "4"), labels = c("Correct model", "Dispersion too large", "Dispersion too small", "Zero inflation")),
         method = factor(x = as.character(method), levels = c("sceptre", "negative_binomial", "scMAGeCK"), labels = c("SCEPTRE", "Fixed dispersion NB", "scMAGeCK")))

p <- qq_data %>%
  ggplot(aes(x = expected, col = method, y = pvalue, ymin = clower, ymax = cupper)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 1, alpha = 0.5) +
  xlab("Expected null p-value") +
  ylab("Observed p-value") +
  scale_colour_manual(values = setNames(plot_colors[c("sceptre", "fixed_dispersion_nb", "scMAGeCK")], NULL), name = "Method") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  facet_wrap(.~facet_label, nrow = 2) +
  theme_bw() + theme(
    panel.spacing.x = unit(1.25, "lines"),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line())

ggsave(plot = p, filename = paste0(offsite_dir, "/figures/simulation_qqplots.pdf"), device = "pdf", scale = 1, width = 5, height = 3.5)