args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
require(katsevich2020)
require(ggrepel)
require(cowplot)
library(dplyr)
library(fst)
library(tidyr)
source(paste0(code_dir, "/sceptre_paper/plotting/load_data_for_plotting.R"))
fig5_dir <- paste0(manuscript_figure_dir, "/Figure5")

