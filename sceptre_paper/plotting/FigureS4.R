# Reproduce Figure S3 from Katsevich, Barry, and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
require(katsevich2020)
require(ggrepel)
require(cowplot)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
figS4_dir <- paste0(manuscript_figure_dir, "/FigureS4")

resampling_results <- resampling_results_xie_cis
original_results <- ss_xie_cis %>% select('gene_id', 'gRNA_id', 'ss.down', 'reject.down') %>% dplyr::rename(rejected = reject.down)
original_results <- left_join(original_results, resampling_results[, c('gene_id', 'gRNA_id', 'gene_short_name', 'chr', 'strand', 'target_gene.start',
                                                   'target_gene.stop', 'TSS', 'target_site.start', 'target_site.stop',
                                                   'target_site.mid')], 
          by = c("gene_id", "gRNA_id"))

# subfigure a (distances)
df = original_results %>% 
  dplyr::rename(rejected_old = rejected, old_pvalue = ss.down) %>%
  #filter(site_type == "cis") %>%
  left_join(resampling_results %>%
              dplyr::select(gene_id, gRNA_id, rejected) %>%
              dplyr::rename(rejected_new = rejected),
            by = c("gene_id", "gRNA_id")) %>%
  mutate(strand = ifelse(target_gene.start == TSS, "+", "-"),
         enhancer_location = 0.5*(target_site.start + target_site.stop),
         TSS_dist = ifelse(strand == "+", enhancer_location - TSS, TSS - enhancer_location)) %>%
  mutate(rejected_new_unique = rejected_new & !rejected_old, 
         rejected_old_unique = rejected_old & !rejected_new)

df1 = df %>% filter(rejected_old_unique) %>% mutate(group = "old_rejections")
df2 = df %>% filter(rejected_new_unique) %>% mutate(group = "new_rejections")

dist_tab <- rbind(
  df1 %>%
    filter(TSS_dist <= 0) %>%
    summarise(`Mean distance (kb)` = mean(TSS_dist)/1000,
              `Median distance (kb)` = median(TSS_dist)/1000) %>%
    mutate(method = "Original"),
  df2 %>%
    filter(TSS_dist <= 0) %>%
    summarise(`Mean distance (kb)` = mean(TSS_dist)/1000, `Median distance (kb)` = median(TSS_dist)/1000) %>%
    mutate(method = "SCEPTRE")
) %>%
  gather(metric, measure, -method) %>%
  spread(method, measure) # %>% kable(format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",digits = 1, col.names = c("", "Original", "SCEPTRE"))
#metric  Original   SCEPTRE
#1   Mean distance (kb) -390.1217 -173.2765
#2 Median distance (kb) -476.2870 -106.9480

p_a <- rbind(df1, df2) %>%
  filter(TSS_dist <= 0) %>%
  mutate(group = factor(group,
                        levels = c("new_rejections", "old_rejections"),
                        labels = c("SCEPTRE", "Virtual FACS"))) %>%
  mutate(TSS_dist = TSS_dist/1000) %>%
  filter(TSS_dist >= -300) %>% arrange(group) %>% 
  ggplot(aes(x = TSS_dist, group = desc(group), fill = group)) +
  geom_histogram(position = position_identity(), alpha = 0.75, bins = 10, aes(y = ..density..)) +
  scale_fill_manual(values = c(plot_colors[["sceptre"]], plot_colors[["hypergeometric"]])) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), limits = c(0, 0.01)) +
  xlab("Enhancer to TSS distance (kb)") + ylab("Frequency") + ggtitle("Gene-enhancer distances") +
  theme_bw() + theme(panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     legend.title = element_blank(),
                     legend.background = element_blank(),
                     legend.position = c(0.2, 0.8),
                     plot.title = element_text(hjust = 0.5))

# subfigure b (TADs)
old_scores = rejected_pairs_HIC_xie %>% filter(rejected_vf & !rejected_sceptre, 
                                                !is.na(score_rank)) %>% pull(score_rank)
new_scores = rejected_pairs_HIC_xie %>% filter(rejected_sceptre & !rejected_vf, 
                                                !is.na(score_rank)) %>% pull(score_rank)

sceptre_percent = rejected_pairs_HIC_xie %>% 
  filter(!rejected_vf & rejected_sceptre) %>%
  summarise(total = dplyr::n(), 
            same_tad = sum(!is.na(TAD_left)), 
            same_tad_prop = sum(!is.na(TAD_left))/dplyr::n())

xie_percent = rejected_pairs_HIC_xie %>% 
  filter(rejected_vf & !rejected_sceptre) %>%
  summarise(total = dplyr::n(), 
            same_tad = sum(!is.na(TAD_left)), 
            same_tad_prop = sum(!is.na(TAD_left))/dplyr::n())


x = seq(0, 1 + 1e-10, length.out = 1000)
old_ecdf = ecdf(old_scores)(x)
new_ecdf = ecdf(new_scores)(x)
p_b <- tibble(x, old_ecdf, new_ecdf) %>%
  gather(method, ecdf, -x) %>%
  mutate(method = factor(method, levels = c("new_ecdf", "old_ecdf"), labels = c("SCEPTRE", "Virtual FACS"))) %>%
  ggplot(aes(x = x, y = ecdf, group = method, colour = method)) +
  geom_line(lwd = 1.2) + geom_abline(slope = 1, linetype = "dashed") +
  scale_colour_manual(values = c( plot_colors[["sceptre"]], plot_colors[["hypergeometric"]]), name = "Loci pairs") + ggtitle("Gene-enhancer HIC interactions")  +
  xlab("Rank of pair by TAD interaction frequency") +
  ylab("Cumulative fraction of pairs") +
  theme_bw() + theme(legend.position = c(0.2, 0.72),
                     legend.title = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     axis.line = element_line(),
                     plot.title = element_text(hjust = 0.5))

# subfigure c (Chip-seq)
my_order <- TF_enrichments_xie %>% filter(method == "SCEPTRE unique Virtual FACS") %>% pull(enrichment) %>% order()
ordered_labs <- (TF_enrichments_xie %>% filter(method == "SCEPTRE unique Virtual FACS") %>% pull(TF))[my_order]
p_c <- TF_enrichments_xie %>% arrange(desc(method)) %>% 
  filter(method %in% c("Virtual FACS unique", "SCEPTRE unique Virtual FACS")) %>%
  mutate(method = factor(method, levels = c("Virtual FACS unique", "SCEPTRE unique Virtual FACS"), 
                         labels = c("Virtual FACS", "SCEPTRE")),
         TF = factor(TF,
                     levels = ordered_labs,
                     labels = ordered_labs)) %>%
  ggplot(aes(x = TF, y = enrichment, fill = method)) +
  geom_col(position = "dodge", width = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  xlab("ChIP-seq target") + ylab("Enrichment (odds ratio)") + ggtitle("Enhancer ChIP-seq enrichment") +
  scale_fill_manual(values = c(plot_colors[["hypergeometric"]], plot_colors[["sceptre"]])) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + coord_flip() + theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, vjust = 0.2, hjust=0))

# combine plots: option 1
combined_xie = plot_grid(p_a, p_b, p_c, ncol = 1, 
                          labels = c("a", "b", "c"), align = "vh", hjust = -4.5)

ggsave(filename = paste0(figS4_dir, "/FigureS4.pdf"), plot = combined_xie, device = "pdf", scale = 1, width = 4.0, height = 7)

