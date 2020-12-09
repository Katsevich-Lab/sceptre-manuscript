# Download data
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_xie/paths_to_dirs.R"))
suppressPackageStartupMessages(library(R.utils))

################################
# 1. Mosaic-seq single-cell data
################################
# Set the source and the destination; perform download and untar
dest <- paste0(raw_data_dir, "/GSE129837_RAW.tar")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129837&format=file", destfile = dest)
untar(dest, exdir = raw_data_dir)
file.remove(dest)

############################
# 2. Bulk RNA-seq validation
############################
dest <- paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.ARL15_enhancer.txt.gz")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129825&format=file&file=GSE129825%5FLibraries%2EFeatureCounts%2EARL15%5Fenhancer%2Etxt%2Egz", destfile = dest)

dest <- paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.MYB_enhancer.txt.gz")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129825&format=file&file=GSE129825%5FLibraries%2EFeatureCounts%2EMYB%5Fenhancer%2Etxt%2Egz", destfile = dest)

# Unzip all the txt.gz files
all_file_names <- list.files(raw_data_dir)
to_unzip <- paste0(raw_data_dir, "/", grep(pattern = '*.gz', x = all_file_names, value = TRUE))
for (file in to_unzip) {
  if (file.exists(file)) {
    gunzip(file)
  }
}

############################
# 3. Spreadsheets and tables
############################
dest <- paste0(raw_data_dir, "/all_oligos.xlsx")
download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124719313956-mmc2.xlsx", destfile = dest)

dest <- paste0(raw_data_dir, "/enh_targets.xlsx")
download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124719313956-mmc4.xlsx", destfile = dest)

dest <- paste0(raw_data_dir, "/bulk_rna_info.xlsx")
download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124719313956-mmc3.xlsx", destfile = dest)

# We put the Genes.xls file in Google drive because automatic download did not work. Download from this link. (source: Human protein-coding genes and gene feature statistics in 2019 by Piovesan et al in BMC Research Notes). 
dest <- paste0(raw_data_dir, "/Genes.xlsx")
download.file(url = "https://doc-0c-48-docs.googleusercontent.com/docs/securesc/njbs04o6ok5oautf8aa0vr10vsh9uq1m/d28vr6jsngqrljapm95kehego46jou0g/1607540850000/15855371786580792779/15855371786580792779/1_6qH8T_m_eX50Ky9y80vQAQwAAKJOP5b?e=download&authuser=0&nonce=dpd6m4u8mbmck&user=15855371786580792779&hash=m61bmqvoj8m9", destfile = dest)

###############################
# 4. Results from other authors
###############################
dest <- paste0(raw_data_dir, "/hypergeometric_pvals_arl15_down.mat")
download.file(url = "https://github.com/russellxie/Global-analysis-K562-enhancers/blob/master/Notebooks/Data/Hypergeometric_pvals/chr5-54325645-54326045-down_log-pval.mat?raw=true", destfile = dest)

dest <- paste0(raw_data_dir, "/hypergeometric_pvals_myb3_down.mat")
download.file(url = "https://github.com/russellxie/Global-analysis-K562-enhancers/blob/master/Notebooks/Data/Hypergeometric_pvals/chr6-135323137-135323537-down_log-pval.mat?raw=true", destfile = dest)
