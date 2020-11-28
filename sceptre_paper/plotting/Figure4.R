# Reproduce Figure 4 from Katsevich, Barry, and Roeder (2020).
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper" else args[1]
require(katsevich2020)
require(ggrepel)
require(cowplot)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
fig4_dir <- paste0(manuscript_figure_dir, "/Figure4")

# subfigure a: old vs. new rejections
annotations <- tibble(gene_short_name = c("B3GNT2", "PTPN1", "TOP1",  "AGFG1", "EIF1"),
                     eQTL = c(2.7e-26, NA, NA, 5.2e-8,NA),
                     eRNA = c(NA, 2e-18, 6.6e-5,NA, 1.2e-6))

promising_new_rejections <- resampling_results_gasp %>%
  filter(gene_short_name %in% c("B3GNT2", "PTPN1", "TOP1",  "AGFG1", "EIF1"), rejected) %>%
  inner_join(original_results_gasp %>% filter(pvalue.empirical.adjusted > 0.1) %>%
               select(gene_id, grna_group, pvalue.empirical), by = c("gene_id", "grna_group")) %>%
  left_join(annotations, by = "gene_short_name") %>%
  arrange(desc(quality_rank_grna), p_value) %>%
  mutate(pair_number = row_number()) %>%
  select(pair_number, pair_id, gene_short_name, target_site, p_value,
         pvalue.empirical, eQTL, eRNA, quality_rank_grna)


thresh_new <- resampling_results_gasp %>% filter(rejected) %>% summarise(max(p_value)) %>% pull()
thresh_old <- original_results_gasp %>% filter(rejected) %>% summarise(max(pvalue.empirical)) %>% pull()
df_a <- resampling_results_gasp %>%
  rename(new_pvalue = p_value, new_rejected = rejected) %>%
  select(pair_id, gene_id, grna_group, new_pvalue, site_type, quality_rank_grna, new_rejected) %>%
  left_join(original_results_gasp %>%
              rename(old_pvalue = pvalue.empirical, old_rejected = rejected) %>%
              select(gene_id, grna_group, old_pvalue, old_rejected, beta, outlier_gene), by = c("gene_id", "grna_group")) %>%
  filter(site_type == "DHS", quality_rank_grna == "top_two") %>%
  select(-c("site_type", "quality_rank_grna")) %>%
  mutate(reason_not_rejected = ifelse(outlier_gene,
                                      "Outlier gene",
                                      ifelse(beta > 0,
                                             "Positive effect",
                                             ifelse(!old_rejected & old_pvalue <= thresh_old, "Low confidence", "none"))),
         reason_not_rejected = factor(reason_not_rejected, levels = c("none", "Positive effect", "Outlier gene", "Low confidence")),
         rejected_by = ifelse(old_rejected,
                              ifelse(new_rejected, "Both methods", "Gasperini only"),
                              ifelse(new_rejected, "SCEPTRE only", "Neither method"))) %>%
  left_join(promising_new_rejections %>%
              select(pair_number, pair_id) %>%
              mutate(pair_number = as.character(pair_number)), by = "pair_id") %>%
  mutate(pair_number = ifelse(is.na(pair_number), "", pair_number)) %>%
  arrange(old_rejected, new_rejected) %>%
  filter(old_pvalue < 1e-2 |
           new_pvalue < 1e-2 |
           reason_not_rejected != "none" |
           row_number() %% 10 == 0)
n_rejected <- table(df_a$rejected_by)
n_rejected_annotated <- ifelse(names(n_rejected) == "Neither method", "Neither method",  paste0(names(n_rejected), " (", n_rejected,")"))
matches <- match(x = df_a$rejected_by, table = names(n_rejected))
rejected_by_annotated <- n_rejected_annotated[matches]
my_labs <- n_rejected_annotated[match(c("Neither method", "Both methods", "SCEPTRE only", "Gasperini only"), names(n_rejected))]
df_a$rejected_by_annotated <- factor(x = rejected_by_annotated, levels = my_labs, labels = my_labs)
df_a$outlier_gene <- as.logical(df_a$outlier_gene)
filter(df_a, outlier_gene, new_rejected) %>% pull(gene_id) %>% unique() %>% length()
filter(df_a, outlier_gene, new_rejected) %>% nrow()

p_a <- ggplot(data = arrange(df_a, rejected_by_annotated), aes(x = old_pvalue, y = new_pvalue,
             color = rejected_by_annotated, label = pair_number, shape = outlier_gene)) + geom_point(alpha = 0.9) +
  geom_text_repel(colour = "black", force = 0.8, point.padding = 0.1, box.padding = 0.4, size = 4, min.segment.length = 0) +
  geom_vline(xintercept = thresh_old, linetype = "dashed") + geom_hline(yintercept = thresh_new, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("grey60", "darkorchid1",  plot_colors[["sceptre"]], plot_colors[["gasperini_nb"]]), name = "Discovered by") +
  scale_shape_manual(values = c("circle", "triangle"), labels = c("Non-outlier", "Outlier"),name = "Gasperini analysis") +
  scale_x_log10() + scale_y_log10() +
  xlab("Gasperini empirical p-value for gene-gRNA pair") +
  guides(color = guide_legend(order = 1)) +
  guides(shape = guide_legend(order = 2)) +
  ylab("SCEPTRE p-value for gene-gRNA pair") + ggtitle("SCEPTRE vs. Gasperini p-values") +
  theme_bw() + theme(
    legend.position = c(0.77, 0.22),
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.background = element_rect(color = "transparent", fill = "white"),
    legend.margin = margin(-0.25,0,0,0, unit="cm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(hjust = 0.5))

# subfigure b (distances)
df = original_results_gasp %>%
  rename(old_rejected = rejected, old_pvalue = pvalue.empirical) %>%
  filter(site_type == "DHS", quality_rank_grna == "top_two") %>%
  left_join(resampling_results_gasp %>%
              select(gene_id, grna_group, rejected) %>%
              rename(new_rejected = rejected),
            by = c("gene_id", "grna_group")) %>%
  mutate(strand = ifelse(target_gene.start == TSS, "+", "-"),
         enhancer_location = 0.5*(target_site.start + target_site.stop),
         TSS_dist = ifelse(strand == "+", enhancer_location - TSS, TSS - enhancer_location))

df1 = df %>% filter(old_rejected) %>% mutate(group = "old_rejections")
df2 = df %>% filter(new_rejected) %>% mutate(group = "new_rejections")

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

p_b <- rbind(df1, df2) %>%
  filter(TSS_dist <= 0) %>%
  mutate(group = factor(group,
                        levels = c("new_rejections", "old_rejections"),
                        labels = c("SCEPTRE", "Original"))) %>%
  mutate(TSS_dist = TSS_dist/1000) %>%
  filter(TSS_dist >= -300) %>%
  ggplot(aes(x = TSS_dist, group = group, fill = group)) +
  geom_histogram(position = position_identity(), alpha = 0.75, bins = 14) +
  scale_fill_manual(values = c(plot_colors[["sceptre"]], plot_colors[["gasperini_nb"]])) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  xlab("Enhancer to TSS distance (kb)") + ylab("Count") + ggtitle("Gene-enhancer distances") +
  theme_bw() + theme(panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     legend.title = element_blank(),
                     legend.background = element_blank(),
                     legend.position = c(0.15, 0.72),
                     plot.title = element_text(hjust = 0.5))

# subfigure c (TADs)
old_scores = rejected_pairs_HIC %>% filter(rejected_old, !is.na(score_rank)) %>% pull(score_rank)
new_scores = rejected_pairs_HIC %>% filter(rejected_new, !is.na(score_rank)) %>% pull(score_rank)

sceptre_percent <- length(new_scores)/(df_a$new_rejected %>% sum)
gasp_percent <- length(old_scores)/(df_a$old_rejected %>% sum)

x = seq(0, 1 + 1e-10, length.out = 1000)
old_ecdf = ecdf(old_scores)(x)
new_ecdf = ecdf(new_scores)(x)
p_c <- tibble(x, old_ecdf, new_ecdf) %>%
  gather(method, ecdf, -x) %>%
  mutate(method = factor(method, levels = c("new_ecdf", "old_ecdf"), labels = c("SCEPTRE", "Original"))) %>%
  ggplot(aes(x = x, y = ecdf, group = method, colour = method)) +
  geom_line(lwd = 1.2) + geom_abline(slope = 1, linetype = "dashed") +
  scale_colour_manual(values = c( plot_colors[["sceptre"]], plot_colors[["gasperini_nb"]]), name = "Loci pairs") + ggtitle("Gene-enhancer TAD pairings")  +
  xlab("Rank of pair by TAD interaction frequency") +
  ylab("Cumulative fraction of pairs") +
  theme_bw() + theme(legend.position = c(0.15, 0.72),
                     legend.title = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     axis.line = element_line(),
                     plot.title = element_text(hjust = 0.5))

# subfigure d (Chip-seq)
my_order <- TF_enrichments %>% filter(method == "Proposed") %>% pull(enrichment) %>% order()
ordered_labs <- (TF_enrichments %>% filter(method == "Proposed") %>% pull(TF))[my_order]

p_d <- TF_enrichments %>% arrange(desc(method)) %>%
  mutate(method = factor(method, levels = c("Original", "Proposed"), labels = c("Original", "SCEPTRE")),
         TF = factor(TF,
                     levels = ordered_labs,
                     labels = ordered_labs)) %>%
  ggplot(aes(x = TF, y = enrichment, fill = method)) +
  geom_col(position = "dodge") +
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

# subfigure e (GTEx eQTL FANTOM)
# make a table of the selected pairs with title
p_e <- ggplot() + theme_bw() + ggtitle("eQTL and eRNA validation") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

# combine plots: option 1
left_col <- plot_grid(p_a, p_d, align = "v", labels = c("a", "e"), ncol = 1, rel_heights = c(2,1))
right_col <- plot_grid(p_e, p_b, p_c, align = "vh", labels = c("b", "c", "d"), ncol = 1)
final_plot <- plot_grid(left_col, right_col, align = "h", ncol = 2, rel_widths = c(1.1,1))
ggsave(filename = paste0(fig4_dir, "/subfigures_a_thru_e.pdf"), plot = final_plot, device = "pdf", scale = 1, width = 8, height = 8)
