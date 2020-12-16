require(tidyverse)
packages <- c("Seurat", "R.utils", "fst", "sn", "MASS", "VGAM", "tidyverse", "bigstatsr", "monocle", "sceptre", "katsevich2020", "furrr", "Seurat", "R.matlab", "R.utils", "fst", "sn", "MASS", "VGAM", "tidyverse", "bigstatsr", "openxlsx", "rhdf5", "sceptre", "katsevich2020", "ggpubr", "cowplot")
packages <- packages %>% unique %>% sort
for (pack in packages) {
  library(pack, character.only = TRUE)
}
map(set_names(packages, packages), function(i) packageVersion(i) %>% as.character) %>% unlist()
