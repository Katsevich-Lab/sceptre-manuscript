#' Reverse log transformation
#'
#' @param base base of the log (default e)
#' @export
revlog_trans <- function(base = exp(1)) {
  ## Define the desired transformation.
  trans <- function(x) {
    -log(x, base)
  }
  ## Define the reverse of the desired transformation
  inv <- function(x){
    base^(-x)
  }
  ## Creates the transformation
  trans_new(paste("revlog-", base, sep = ""),
            trans, ## The transformation function (can be defined using anonymous functions)
            inv,  ## The reverse of the transformation
            log_breaks(base = base), ## default way to define the scale breaks
            domain = c(1e-100, Inf) ## The domain over which the transformation is valued
  )
}


#' Plot skew-t given gene and gRNA.
#'
#' @param gene_id id of gene
#' @param gRNA_id id of gRNA
#' @param sceptre_function_args path to the sceptre_function_args.R file. Defaults to the location on Tim's machine.
#' @param offsite_dirpath to the offsite dir. Defaults to the location on Tim's machine.
#'
#' @return a ggplot containing the plot
#' @export
#'
#' @examples
#' gene_id <- "ENSG00000008256"
#' gRNA_id <- "ACTB_TSS"
#' plot_skew_t_gene_gRNA(gene_id, gRNA_id)
plot_skew_t_gene_gRNA <- function(gene_id, gRNA_id, sceptre_function_args = "/Users/timbarry/Box/SCEPTRE/sceptre_paper/sceptre_paper/analysis_drivers/analysis_drivers_gasp/sceptre_function_args.R", offsite_dir = "/Volumes/tims_new_drive/research/sceptre_files") {
  info_pack <- get_all_data_for_gene_gRNA(gene_id, gRNA_id, sceptre_function_args, offsite_dir)
  sceptre_res <- run_sceptre_using_precomp(expressions = info_pack$expressions, gRNA_indicators = info_pack$gRNA_indicators, gRNA_precomp = info_pack$gRNA_precomp, gene_precomp_size = info_pack$gene_precomp_size, gene_precomp_offsets = info_pack$gene_precomp_offsets, B = 500,seed = 1234, reduced_output = FALSE)
  p <- plot_skew_t(resampled_zvalues = sceptre_res$resampled_z_values, original_zvalue = sceptre_res$z_value, dp = sceptre_res$skew_t_mle)
  return(list(plot = p, p_val = sceptre_res$p_value, skew_t_fit_success = sceptre_res$skew_t_fit_success))
}


#' Plot skew-t
#'
#' @param resampled_zvalues the set of resampled z-values
#' @param original_zvalue the z value of the original negative binomial fit
#' @param dpv the skew-t fit MLE
#'
#' @return a ggplot object containing the plot
plot_skew_t <- function(resampled_zvalues, original_zvalue, dp) {
  z = seq(-4, 4, length.out = 1000)
  df_curves = tibble(z = z, fitted = dst(x = z, dp = dp), gaussian = dnorm(z)) %>%
    gather("curve", "y", fitted, gaussian) %>%
    mutate(curve = factor(curve, levels = c("fitted","gaussian"), labels = c("Conditional\nrandomization","Negative\nbinomial")))
  df_ribbon = df_curves %>%
    filter(z <= original_zvalue, curve == "Conditional\nrandomization") %>%
    mutate(lower = 0, upper = y) %>%
    select(z, lower, upper)
  p <- ggplot() +
    geom_histogram(aes(x = z, y = ..density..),
                   data = tibble(z = resampled_zvalues),
                   boundary = 0, colour = "black", fill = "lightgray", binwidth = 0.5) +
    geom_line(aes(x = z, y = y, group = curve, colour = curve, linetype = curve), data = df_curves) +
    geom_vline(xintercept = original_zvalue, colour = "firebrick3", linetype = "solid") +
    geom_ribbon(aes(x = z, ymin = lower, ymax = upper), fill = "darkorchid2", alpha = 0.5, data = df_ribbon) +
    scale_colour_manual(values = c("darkorchid2", "black"), name = "Null distribution") +
    scale_linetype_manual(values = c("solid", "dashed"), name = "Null distribution") +
    scale_y_continuous(expand = c(0,0)) +
    xlab("Negative binomial z-value") + theme_bw() +
    theme(legend.position = c(0.85,0.8),
          legend.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(hjust = 0.5, size = 11),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  return(p)
}


#' Get all data for gene-gRNA pair on Gasperini data
#'
#' @param gene_id id of the gene
#' @param gRNA_id od of the gRNA
#' @param sceptre_function_args path to the sceptre_function_args.R file. Defaults to the location on Tim's machine.
#' @param offsite_dir path to the offsite dir. Defaults to the location on Tim's machine.
#'
#' @return a list containing expressions, gRNA indicators, gRNA precomp, gene precomp size, gene precomp offsets.
#' @export
#'
#' @examples
#' gene_id <- "ENSG00000008256"
#' gRNA_id <- "ACTB_TSS"
#' get_all_data_for_gene_gRNA(gene_id, gRNA_id)
get_all_data_for_gene_gRNA <- function(gene_id, gRNA_id, sceptre_function_args = "/Users/timbarry/Box/SCEPTRE/sceptre_paper/sceptre_paper/analysis_drivers/analysis_drivers_gasp/sceptre_function_args.R", offsite_dir = "/Volumes/tims_new_drive/research/sceptre_files") {
  source(sceptre_function_args)
  expressions <- cell_gene_expression_matrix[, which(ordered_gene_ids == gene_id)][cell_subset]
  gRNA_indicators <- (read.fst(gRNA_indicator_matrix_fp, columns = gRNA_id) %>% pull() %>% as.integer())[cell_subset]
  gene_sizes <- readRDS(paste0(gene_precomp_dir, "/size_reg_file.rds"))
  gene_precomp_size <- gene_sizes[[gene_id]]
  gRNA_precomp <- paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst") %>% read.fst() %>% filter(id == gRNA_id) %>% pull(pod_id) %>% paste0(gRNA_precomp_dir, "/gRNA_precomp_", ., ".fst") %>% read.fst(columns = gRNA_id) %>% pull()
  gene_precomp_offsets <- paste0(gene_precomp_dir, "/gene_dictionary.fst") %>% read.fst() %>% filter(id == gene_id) %>% pull(pod_id) %>% paste0(gene_precomp_dir, "/gene_offsets_", ., ".fst") %>% read.fst(columns = gene_id) %>% pull()
  return(list(expressions = expressions, gRNA_indicators = gRNA_indicators, gRNA_precomp = gRNA_precomp, gene_precomp_size = gene_precomp_size, gene_precomp_offsets = gene_precomp_offsets))
}
