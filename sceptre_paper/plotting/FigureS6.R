p_thresh <- 1e-8
qq_data <- simulation_results %>%
  group_by(method, dataset_id) %>%
  dplyr:::mutate(r = rank(p_value, ties.method = "first"), expected = ppoints(dplyr::n())[r],
                 clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = dplyr::n()+1-r),
                 cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = dplyr::n()+1-r)) %>%
  ungroup() %>%
  dplyr:::mutate(pvalue = ifelse(p_value < p_thresh, p_thresh, p_value),
                 facet_label = factor(x = as.character(dataset_id), levels = c("2", "1", "3", "4"), labels = c("Correct model", "Dispersion too large", "Dispersion too small", "Zero inflation")),
                 method = factor(x = as.character(method), levels = c("sceptre", "negative_binomial", "scMAGeCK"), labels = c("SCEPTRE", "Negative Binomial", "scMAGeCK")))

qq_data <- qq_data[qq_data$method == 'scMAGeCK', ]
p <- qq_data %>%
  ggplot(aes(x = expected, y = pvalue, ymin = clower, ymax = cupper)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 1, alpha = 0.5, color = plot_colors[c("scMAGeCK")]) +
  xlab("Expected null p-value") +
  ylab("Observed p-value") +
  ggtitle("Simulated negative control pair, scMAGeCK") +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  facet_wrap(.~facet_label, nrow = 1) +
  theme_bw() + theme(
    legend.position = c(0.15, 0.8),
    legend.title = element_blank(),
    panel.spacing.x = unit(1.25, "lines"),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line())

ggsave(filename = paste0(figS6_dir, "/FigureS6.pdf"), plot = p, device = "pdf", scale = 1, width = 7, height = 2)
