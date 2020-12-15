#' Plot skew-t
#'
#' @param resampled_zvalues the set of resampled z-values
#' @param original_zvalue the z value of the original negative binomial fit
#' @param dpv the skew-t fit MLE
#'
#' @return a ggplot object containing the plot
#' @export
plot_skew_t <- function(resampled_zvalues, original_zvalue, dp, interval = c(-4,4)) {
  z = seq(interval[1], interval[2], length.out = 1000)
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
                   boundary = 0, colour = "black", fill = "lightgray", bins = 15) +
    geom_line(aes(x = z, y = y, group = curve, colour = curve, linetype = curve), data = df_curves) +
    geom_vline(xintercept = original_zvalue, colour = "firebrick3", linetype = "solid") +
    geom_ribbon(aes(x = z, ymin = lower, ymax = upper), fill = "darkorchid2", alpha = 0.5, data = df_ribbon) +
    scale_colour_manual(values = c("darkorchid2", "black"), name = "Null distribution") +
    scale_linetype_manual(values = c("solid", "dashed"), name = "Null distribution") +
    scale_y_continuous(expand = c(0,0)) +
    xlab("Negative binomial z-value") + theme_bw() +
    theme(legend.position = c(0.85, 0.8),
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
