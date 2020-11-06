args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
suppressPackageStartupMessages(library(katsevich2020))
suppressPackageStartupMessages(library(scales))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
simulation_param_file <- if (is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/simulations/simulation_parameters.R" else args[2]
source(simulation_param_file)

fs <- list.files(simulation_result_dir, full.names = TRUE)
res <- fs %>% map(readRDS) %>% reduce(rbind)

sc1 <- filter(res, theta_size == "theta_small" & method == "sceptre")
sc2 <- filter(res, theta_size == "theta_correct" & method == "sceptre")

ci <- 0.95
truncate_thresh <- 1e-9
qq_data <- res %>%
  group_by(method, theta_size) %>%
  mutate(r = rank(p_value), expected = ppoints(n())[r],
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>%
  mutate(pvalue = ifelse(p_value < truncate_thresh, truncate_thresh, p_value))

p <- qq_data %>%
  ggplot(aes(x = expected, y = pvalue, col = method, ymin = clower, ymax = cupper)) +
  geom_point(aes(color = method), size = 1, alpha = 0.5) +
  geom_ribbon(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Expected null p-value") +
  ylab("Observed p-value") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) + facet_grid(.~theta_size)
