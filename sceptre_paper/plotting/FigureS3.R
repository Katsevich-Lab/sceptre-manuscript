# Reproduce Figure S3 from Katsevich, Barry, and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
require(katsevich2020)
require(ggrepel)
require(cowplot)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
figS3_dir <- paste0(manuscript_figure_dir, "/FigureS3")

# subfigure a (distances)
df = original_results_gasp %>% 
  dplyr::rename(rejected_old = rejected, old_pvalue = pvalue.empirical) %>%
  filter(site_type == "DHS", quality_rank_grna == "top_two") %>%
  left_join(resampling_results_gasp %>%
              dplyr::select(gene_id, grna_group, rejected) %>%
              dplyr::rename(rejected_new = rejected),
            by = c("gene_id", "grna_group")) %>%
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

p_a <- rbind(df1, df2) %>%
  filter(TSS_dist <= 0) %>%
  mutate(group = factor(group,
                        levels = c("new_rejections", "old_rejections"),
                        labels = c("SCEPTRE", "Original"))) %>%
  mutate(TSS_dist = TSS_dist/1000) %>%
  filter(TSS_dist >= -300) %>% arrange(group) %>% 
  ggplot(aes(x = TSS_dist, group = desc(group), fill = group)) +
  geom_histogram(position = position_identity(), alpha = 0.75, bins = 10, aes(y = ..density..)) +
  scale_fill_manual(values = c(plot_colors[["sceptre"]], plot_colors[["gasperini_nb"]])) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  xlab("Enhancer to TSS distance (kb)") + ylab("Frequency") + ggtitle("Gene-enhancer distances") +
  theme_bw() + theme(panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     legend.title = element_blank(),
                     legend.background = element_blank(),
                     legend.position = c(0.17, 0.72),
                     plot.title = element_text(hjust = 0.5))

# subfigure b (TADs)
old_scores = rejected_pairs_HIC %>% filter(rejected_old & !rejected_new, 
                                           !is.na(score_rank)) %>% pull(score_rank)
new_scores = rejected_pairs_HIC %>% filter(rejected_new & !rejected_old, 
                                           !is.na(score_rank)) %>% pull(score_rank)

sceptre_percent = rejected_pairs_HIC %>% 
  filter(!rejected_old & rejected_new) %>%
  summarise(total = dplyr::n(), 
            same_tad = sum(!is.na(TAD_left)), 
            same_tad_prop = sum(!is.na(TAD_left))/dplyr::n())

gasp_percent = rejected_pairs_HIC %>% 
  filter(rejected_old & !rejected_new) %>%
  summarise(total = dplyr::n(), 
            same_tad = sum(!is.na(TAD_left)), 
            same_tad_prop = sum(!is.na(TAD_left))/dplyr::n())

x = seq(0, 1 + 1e-10, length.out = 1000)
old_ecdf = ecdf(old_scores)(x)
new_ecdf = ecdf(new_scores)(x)
p_b <- tibble(x, old_ecdf, new_ecdf) %>%
  gather(method, ecdf, -x) %>%
  mutate(method = factor(method, levels = c("new_ecdf", "old_ecdf"), labels = c("SCEPTRE", "Original"))) %>%
  ggplot(aes(x = x, y = ecdf, group = method, colour = method)) +
  geom_line(lwd = 1.2) + geom_abline(slope = 1, linetype = "dashed") +
  scale_colour_manual(values = c( plot_colors[["sceptre"]], plot_colors[["gasperini_nb"]]), name = "Loci pairs") + ggtitle("Gene-enhancer HIC interactions")  +
  xlab("Rank of pair by TAD interaction frequency") +
  ylab("Cumulative fraction of pairs") +
  theme_bw() + theme(legend.position = c(0.17, 0.72),
                     legend.title = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     axis.line = element_line(),
                     plot.title = element_text(hjust = 0.5))

# subfigure c (Chip-seq)
my_order <- TF_enrichments %>% filter(method == "Proposed Unique") %>% pull(enrichment) %>% order()
ordered_labs <- (TF_enrichments %>% filter(method == "Proposed Unique") %>% pull(TF))[my_order]
p_c <- TF_enrichments %>% arrange(desc(method)) %>%
  mutate(method = factor(method, levels = c("Original Unique", "Proposed Unique"), labels = c("Original", "SCEPTRE")),
         TF = factor(TF,
                     levels = ordered_labs,
                     labels = ordered_labs)) %>%
  ggplot(aes(x = TF, y = enrichment, fill = method)) +
  geom_col(position = "dodge", width = 1.25) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  xlab("ChIP-seq target") + ylab("Enrichment (odds ratio)") + ggtitle("Enhancer ChIP-seq enrichment") +
  scale_fill_manual(values = c(plot_colors[["gasperini_nb"]], plot_colors[["sceptre"]])) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + coord_flip() + theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, vjust = 0.2, hjust=0))

# combine plots: option 1
combined = plot_grid(p_a, p_b, p_c, ncol = 1, 
                     labels = c("a", "b", "c"), align = "vh", hjust = -4.5)
ggsave(filename = paste0(figS3_dir, "/subfigures_a_thru_c.pdf"), plot = combined, device = "pdf", scale = 1, width = 4.0, height = 7)
