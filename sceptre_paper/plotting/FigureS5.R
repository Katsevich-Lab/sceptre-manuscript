# Reproduce Figure S4 from Katsevich and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
require(katsevich2020)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
figS4_dir <- paste0(manuscript_figure_dir, "/FigureS4")

# plot the fraction paired in each quintile
p_a = paired_fractions_gasp %>% 
  gather(method, paired_fraction, rejected_monocle, rejected_sceptre, rejected_improve) %>% 
  mutate(method = factor(method, levels = c("rejected_sceptre", "rejected_monocle", "rejected_improve"), 
                         labels = c("SCEPTRE", "Monocle NB", "Improved NB"))) %>%
  ggplot(aes(x = factor(quintile), y = paired_fraction, fill = method)) + 
  xlab("ChIP-seq quintiles of candidate enhancers") + ylab("Proportion enhancers paired with gene") +
  geom_col(position = "dodge") + scale_fill_manual(values = c(plot_colors[["sceptre"]], plot_colors[["gasperini_nb"]], plot_colors[['hf_nb']])) + 
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(TF ~ ., nrow = 2) + theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.09, 0.88), 
        legend.background = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line()
  )
#plot(p)

p_b = paired_fractions_xie %>% 
  gather(method, paired_fraction, rejected_vf, rejected_monocle, rejected_sceptre, rejected_nb) %>% 
  mutate(method = factor(method, levels = c("rejected_sceptre", "rejected_vf", "rejected_monocle", "rejected_nb"), 
                         labels = c("SCEPTRE", "Virtual FACS", "Monocle NB", "Improved NB"))) %>%
  ggplot(aes(x = factor(quintile), y = paired_fraction, fill = method)) + 
  xlab("ChIP-seq quintiles of candidate enhancers") + ylab("Proportion enhancers paired with gene") +
  geom_col(position = "dodge") + scale_fill_manual(values = c(plot_colors[["sceptre"]], plot_colors[['hypergeometric']], 
                                                              plot_colors[["gasperini_nb"]], plot_colors[['hf_nb']])) + 
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(TF ~ ., nrow = 2) + theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.17, 0.93), 
        legend.background = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line()
  ) + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

p = plot_grid(p_a, p_b, labels = c('a', 'b'), nrow = 2)
ggsave(plot = p, filename = sprintf("%s/FigureS4.pdf", figS4_dir), width = 7, height = 8)
