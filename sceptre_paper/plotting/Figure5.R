df = original_results %>% 
  rename(old_rejected = rejected, old_pvalue = pvalue.empirical) %>%
  filter(site_type == "DHS", quality_rank_grna == "top_two") %>%
  left_join(resampling_results %>% 
              filter(method == "conditional_randomization") %>%
              select(gene_id, grna_group, rejected) %>% 
              rename(new_rejected = rejected), 
            by = c("gene_id", "grna_group")) %>%
  mutate(strand = ifelse(target_gene.start == TSS, "+", "-"),
         enhancer_location = 0.5*(target_site.start + target_site.stop),
         TSS_dist = ifelse(strand == "+", enhancer_location - TSS, TSS - enhancer_location))

df1 = df %>% filter(old_rejected) %>% mutate(group = "old_rejections")
df2 = df %>% filter(new_rejected) %>% mutate(group = "new_rejections")

p = rbind(df1, df2) %>%
  filter(TSS_dist <= 0) %>%
  mutate(group = factor(group, 
                        levels = c("new_rejections", "old_rejections"),
                        labels = c("SCEPTRE", "Original"))) %>%
  mutate(TSS_dist = TSS_dist/1000) %>%
  ggplot(aes(x = TSS_dist, group = group, fill = group)) + 
  geom_histogram(position = position_identity(), alpha = 0.5) + 
  scale_fill_manual(values = c("blue", "red")) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  xlab("Enhancer to TSS distance (kb)") +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     legend.title = element_blank(),
                     legend.position = c(0.2,0.8)) 
plot(p)
ggsave(filename = sprintf("%s/Figure5/TSS_dist_histogram.pdf", figures_dir), plot = p, device = "pdf",
       width = 4, height = 3)

rbind(
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
  spread(method, measure) %>% 
  kable(format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",digits = 1,
        col.names = c("", "Original", "SCEPTRE")) %>%
  save_kable(sprintf("~/Desktop/TSS_dist_summary.pdf"))

old_scores = rejected_pairs_HIC %>% filter(rejected_old, !is.na(score_rank)) %>% pull(score_rank)
new_scores = rejected_pairs_HIC %>% filter(rejected_new, !is.na(score_rank)) %>% pull(score_rank)
x = seq(0, 1 + 1e-10, length.out = 1000)
old_ecdf = ecdf(old_scores)(x)
new_ecdf = ecdf(new_scores)(x)
p = tibble(x, old_ecdf, new_ecdf) %>% 
  gather(method, ecdf, -x) %>% 
  mutate(method = factor(method, levels = c("new_ecdf", "old_ecdf"), labels = c("SCEPTRE", "Original"))) %>%
  ggplot(aes(x = x, y = ecdf, group = method, colour = method)) + 
  geom_line() + geom_abline(slope = 1, linetype = "dashed") + 
  scale_colour_manual(values = c("blue", "red"), name = "Loci pairs") + 
  xlab("Rank of loci pair by HI-C interaction frequency within TAD") + 
  ylab("Cumulative fraction of loci pairs") +
  theme_bw() + theme(
    legend.position = c(0.175,0.825),
    legend.title = element_blank(),
    panel.grid = element_blank(), 
    panel.border = element_blank(), 
    axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/Figure5/Figure5b.pdf", figures_dir), plot = p, device = "pdf",
       width = 5, height = 3)

p = TF_enrichments %>%
  mutate(method = factor(method, levels = c("Proposed", "Original"), labels =c("SCEPTRE", "Original")),
         TF = factor(TF, 
                     levels = c("H3K27ac", "EP300", "BRD4", "GATA2", "TAL1", "TBL1XR1", "DPF2", "RNF2"),
                     labels = c("H3K27ac", "P300", "BRD4", "GATA2", "TAL1", "TBL1XR1", "DPF2", "RNF2"))) %>%
  ggplot(aes(x = TF, y = enrichment, fill = method)) + 
  geom_col(position = "dodge") + 
  xlab("ChIP-seq target") + ylab("Enrichment (odds ratio)\namong paired enhancers") +
  scale_fill_manual(values = c("blue", "red"), name = "Paired enhancers") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() + theme(
    legend.position = "none",
    # legend.position = c(0.825,0.8),
    panel.grid = element_blank(), 
    panel.border = element_blank(), 
    axis.line = element_line())
plot(p)
ggsave(plot = p, filename = sprintf("%s/Figure5/Figure5c.pdf", figures_dir), width = 8, height = 3)
