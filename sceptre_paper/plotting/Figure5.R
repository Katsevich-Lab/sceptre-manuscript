args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
library(mgcv)
library(scales)
require(katsevich2020)
require(ggrepel)
require(cowplot)
library(dplyr)
library(fst)
library(tidyr)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
fig5_dir <- paste0(manuscript_figure_dir, "/Figure5")
####### subfigure a (Significance Score vs SCEPTRE p-values) #########
resampling_results <- resampling_results_xie_cis
original_results <- ss_xie_cis %>% select('gene_id', 'gRNA_id', 'ss.down', 'reject.down') %>% dplyr::rename(rejected = reject.down)
original_results = cbind(original_results, resampling_results[, c('gene_short_name', 'chr', 'strand', 'target_gene.start', 
                                                                  'target_gene.stop', 'TSS', 'target_site.start', 'target_site.stop',
                                                                  'target_site.mid')])

rejected_by_annotated = rep(NA, nrow(resampling_results))
rejected_by_annotated[resampling_results$rejected*10 + ss_xie_cis$reject.down == 0] = 'Neither method'
rejected_by_annotated[resampling_results$rejected*10 + ss_xie_cis$reject.down == 1] = 'Virtual FACS only'
rejected_by_annotated[resampling_results$rejected*10 + ss_xie_cis$reject.down == 10] = 'SCEPTRE only'
rejected_by_annotated[resampling_results$rejected*10 + ss_xie_cis$reject.down == 11] = 'Both methods'
rejected_by_annotated = factor(rejected_by_annotated, levels = c("Neither method", "Both methods", "SCEPTRE only", "Virtual FACS only"))
table(rejected_by_annotated)
#rejected_by_annotated
#   Neither method      Both methods      SCEPTRE only Virtual FACS only 
#             3367                84                51                51 

ss_thres =sort(ss_xie_cis$ss.down, decreasing = T)[135]
df_s = data.frame(gene_id = resampling_results$gene_id, gRNA_id = resampling_results$gRNA_id, 
                  p_value = resampling_results$p_value, signif.score = ss_xie_cis$ss.down, 
                  rejected_by_annotated = rejected_by_annotated)
p_a = ggplot(data = arrange(df_s, rejected_by_annotated), aes(x = signif.score, y = p_value, color = rejected_by_annotated)) + 
  geom_point(alpha = 0.9) +
  #geom_text_repel(colour = "black", force = 0.8, point.padding = 0.1, box.padding = 0.3, size = 3, min.segment.length = 0, max.overlaps = 50) +
  geom_vline(xintercept = ss_thres, linetype = "dashed") + geom_hline(yintercept = 0.1/26, linetype = "dashed") +
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("grey60", "darkorchid1",  plot_colors[["sceptre"]], plot_colors[["hypergeometric"]]), name = "") +
  #scale_shape_manual(values = c("circle", "triangle"), labels = c("Non-outlier", "Outlier"),name = "Gasperini analysis") +
  scale_x_continuous(trans = pseudo_log_trans(base = 10), breaks = c(-5, 0, 10, 100, 1000), limits = c(-7, 9000)) + 
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  xlab("Virtual FACS Significance Score") +
  guides(color = guide_legend(order = 1)) +
  guides(shape = guide_legend(order = 2)) +
  ylab("SCEPTRE p-value") + ggtitle("SCEPTRE p-values vs. Virtual FACS \n Significance Scores") +
  theme_bw() + theme(
    legend.position = c(0.7, 0.45),
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.background = element_blank(),
    legend.margin = margin(-0.25,0,0,0, unit="cm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(hjust = 0.5))

####### subfigure b (distances) ##########
df = original_results %>%
  dplyr::rename(old_rejected = rejected) %>%
  #filter(site_type == "cis") %>%
  left_join(resampling_results %>%
              select(gene_id, gRNA_id, rejected) %>%
              dplyr::rename(new_rejected = rejected),
            by = c("gene_id", "gRNA_id")) %>%
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
#metric  Original  SCEPTRE
#1   Mean distance (kb) -208.6636 -97.8121
#2 Median distance (kb)  -58.5065 -40.9805

p_b <- rbind(df1, df2) %>%
  filter(TSS_dist <= 0) %>%
  mutate(group = factor(group,
                        levels = c("new_rejections", "old_rejections"),
                        labels = c("SCEPTRE", "Virtual FACS"))) %>%
  mutate(TSS_dist = TSS_dist/1000) %>%
  filter(TSS_dist >= -300) %>%
  ggplot(aes(x = TSS_dist, group = desc(group), fill = group)) +
  geom_histogram(position = position_identity(), alpha = 0.75, bins = 14, aes(y = ..density..)) +
  scale_fill_manual(values = c(plot_colors[["sceptre"]], plot_colors[["hypergeometric"]])) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  xlab("Enhancer to TSS distance (kb)") + ylab("Frequency") + ggtitle("Gene-enhancer distances") +
  theme_bw() + theme(panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(),
                     legend.title = element_blank(),
                     legend.background = element_blank(),
                     legend.position = c(0.2, 0.8),
                     plot.title = element_text(hjust = 0.5))

# subfigure c (TADs)
old_scores = rejected_pairs_HIC_xie %>% filter(rejected_old, !is.na(score_rank)) %>% pull(score_rank)
new_scores = rejected_pairs_HIC_xie %>% filter(rejected_new, !is.na(score_rank)) %>% pull(score_rank)
sceptre_percent = rejected_pairs_HIC_xie %>% 
  filter(rejected_new) %>%
  dplyr::summarise(total = dplyr::n(), 
                   same_tad = sum(!is.na(TAD_left)), 
                   same_tad_prop = sum(!is.na(TAD_left))/dplyr::n())

original_percent = rejected_pairs_HIC_xie %>% 
  filter(rejected_old) %>%
  dplyr::summarise(total = dplyr::n(), 
                   same_tad = sum(!is.na(TAD_left)), 
                   same_tad_prop = sum(!is.na(TAD_left))/dplyr::n())
#> sceptre_percent
# A tibble: 1 x 3
#total same_tad same_tad_prop
#<int>    <int>         <dbl>
#  1   135       97         0.719
#> original_percent
# A tibble: 1 x 3
#total same_tad same_tad_prop
#<int>    <int>         <dbl>
#  1   135       83         0.615

x = seq(0, 1 + 1e-10, length.out = 1000)
old_ecdf = ecdf(old_scores)(x)
new_ecdf = ecdf(new_scores)(x)
p_c <- tibble(x, old_ecdf, new_ecdf) %>%
  gather(method, ecdf, -x) %>%
  mutate(method = factor(method, levels = c("new_ecdf", "old_ecdf"), labels = c("SCEPTRE", "Virtual FACS"))) %>%
  ggplot(aes(x = x, y = ecdf, group = method, colour = method)) +
  geom_line(lwd = 1.2) + geom_abline(slope = 1, linetype = "dashed") +
  scale_colour_manual(values = c( plot_colors[["sceptre"]], plot_colors[["hypergeometric"]]), name = "Loci pairs") + ggtitle("Gene-enhancer HIC interactions")  +
  xlab("Rank of pair by TAD interaction frequency") +
  ylab("Cumulative fraction of pairs") +
  theme_bw() + theme(legend.position = c(0.2, 0.8),
                     legend.title = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     axis.line = element_line(),
                     plot.title = element_text(hjust = 0.5))
# subfigure d (Chip-seq)
my_order <- TF_enrichments_xie %>% filter(method == "Proposed") %>% pull(enrichment) %>% order()
ordered_labs <- (TF_enrichments_xie %>% filter(method == "Proposed") %>% pull(TF))[my_order]

p_d <- TF_enrichments_xie %>% arrange(desc(method)) %>%
  mutate(method = factor(method, levels = c("Original", "Proposed"), labels = c("Virtual FACS", "SCEPTRE")),
         TF = factor(TF,
                     levels = ordered_labs,
                     labels = ordered_labs)) %>%
  #filter(method %in% c('Original', 'SCEPTRE')) %>%
  ggplot(aes(x = TF, y = enrichment, fill = method)) +
  geom_col(position = "dodge", width = 1.25) +
  xlab("ChIP-seq target") + ylab("Enrichment (odds ratio)") + ggtitle("Enhancer ChIP-seq enrichment") +
  scale_fill_manual(values = c(plot_colors[["hypergeometric"]], plot_colors[["sceptre"]])) +
  scale_y_continuous(expand = c(0, 0)) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() + coord_flip() + theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, vjust = 0.2, hjust=0))

combined_full = plot_grid(p_a, p_b, p_c, p_d, align = 'v', labels = c("a", "b", "c", "d"), nrow = 2)
ggsave(filename = paste0(fig5_dir, "/Fig5_full.pdf"), plot = combined_full, device = "pdf", scale = 1, width = 8, height = 6.5)
