####################################################################
# Download raw data from the web
# Note: The GeneHancer database is proprietary and therefore
# must be accessed piecemeal via https://genealacart.genecards.org/.
####################################################################
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_gasp/file_paths_to_dirs.R"))
suppressPackageStartupMessages(library(R.utils))
################################
# Download Gasperini et al. data
################################

# path to store the raw Gasperini data
raw_data_dir_gasp <- raw_data_dir

# URL of data
remote <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120861&format=file&file="

# Gasperini et al results
all_deg_results_filename <- "GSE120861_all_deg_results.at_scale.txt"

# names of genes
genes_filename <- "GSE120861_at_scale_screen.genes.txt"

# names of cells
cells_filename <- "GSE120861_at_scale_screen.cells.txt"

# "reference cells" used by Gasperini et al for computational purposes
reference_cells_filename <- "GSE120861_50k_reference_cells.rds"

# all (gRNA, gene) pairs
gRNAgroup_pair_table_filename <- "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt"

# list of gRNA groups used
gRNA_groups_filename <- "GSE120861_grna_groups.at_scale.txt"

# Monocle Cell Data Set object with all gRNA data
cds_filename <- "GSE120861_at_scale_screen.cds.rds"

# Expression data
expression_filename <- "GSE120861_at_scale_screen.exprs.mtx"

# list of files to download
filenames <- c(all_deg_results_filename,
              genes_filename,
              cells_filename,
              reference_cells_filename,
              cds_filename,
              expression_filename,
              gRNAgroup_pair_table_filename,
              gRNA_groups_filename)

# download files if not already present
for (filename in filenames) {
  if (!file.exists(paste0(raw_data_dir_gasp, "/", filename))) {
    cat(paste0("Downloading ", filename, "\n"))
    source <- paste0(remote, filename, ".gz")
    dest <- paste0(raw_data_dir_gasp, "/", filename, ".gz")
    download.file(source, dest)
    gunzip(paste0(dest))
  }
}

# Download supplementary Table S2 from Cell
supplementary_table_file <- "https://www.cell.com/cms/10.1016/j.cell.2018.11.029/attachment/7319ccb0-a8c0-45f3-8203-26b9159b0102/mmc2.xlsx"
download.file(supplementary_table_file, paste0(raw_data_dir_gasp, "/Gasperini_TableS2.xlsx"))

####################
# Download HI-C data
####################
HIC_dir <- paste0(offsite_dir, "/data/functional/HIC")
HIC_loc <- paste0(HIC_dir, "/GSE63525_K562_Arrowhead_domainlist.txt")
download.file(url = "https://www.dropbox.com/sh/65twliq3qt78ex0/AABpoTxeQHJNbG7TsKELmuEwa/data/raw/HIC/GSE63525_K562_Arrowhead_domainlist.txt?dl=1", destfile = HIC_loc)

remote <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63525&format=file&file="
contact_matrices_dirname <- "GSE63525_K562_intrachromosomal_contact_matrices"
if(!dir.exists(paste0(HIC_dir, "/" , contact_matrices_dirname))) {
  cat(sprintf("Downloading %s...\n", contact_matrices_dirname))
  download.file(sprintf("%s%s.tar.gz", remote, contact_matrices_dirname),
                sprintf("%s/%s.tar.gz", HIC_dir, contact_matrices_dirname))
  gunzip(sprintf("%s/%s.tar.gz", HIC_dir, contact_matrices_dirname))
  untar(tarfile = sprintf("%s/%s.tar", HIC_dir, contact_matrices_dirname), exdir = sprintf("%s/%s", HIC_dir, contact_matrices_dirname))
}

########################
# Download ChIP-seq data (Note: these data were delivered by the Shendure Lab. We put them into Dropbox so that they are permanently downloadable.)
########################
chip_seq_files <- c("BRD4", "DPF2", "EP300", "GATA2", "H3K27ac", "RNF2", "TAL1", "TBL1XR1")
chip_seq_file_location <- "bit.ly/SCEPTRE/raw/ChIP-seq"
# download.file(url = sprintf("%s/BRD4.bed?dl=1", chip_seq_file_location), destfile = "./test")
# download.file(url = urls[i], destfile = dest_locs[i])

urls <- c("https://www.dropbox.com/sh/65twliq3qt78ex0/AABu0Vgsj4Wf-8dijc_mIDssa/data/raw/ChIP-seq/BRD4.bed?dl=1",
          "https://www.dropbox.com/sh/65twliq3qt78ex0/AACLMKyUAsNFudCTHwZpu0eja/data/raw/ChIP-seq/DPF2.bed?dl=1",
          "https://www.dropbox.com/sh/65twliq3qt78ex0/AAAsMLrZLkEf6Q5T0wTSmiLwa/data/raw/ChIP-seq/EP300.bed?dl=1",
          "https://www.dropbox.com/sh/65twliq3qt78ex0/AADtl80uTvunTt1CNUd_4lUra/data/raw/ChIP-seq/GATA2.bed?dl=1",
          "https://www.dropbox.com/sh/65twliq3qt78ex0/AABnS8MAk4-gBljbWG1fpgGJa/data/raw/ChIP-seq/H3K27ac.bed?dl=1",
          "https://www.dropbox.com/sh/65twliq3qt78ex0/AAAGRES3MsckonN8ZAwr_8jYa/data/raw/ChIP-seq/RNF2.bed?dl=1",
          "https://www.dropbox.com/sh/65twliq3qt78ex0/AACIaHcTumi85sUGqmkH78Gsa/data/raw/ChIP-seq/TAL1.bed?dl=1",
          "https://www.dropbox.com/sh/65twliq3qt78ex0/AADVYs4mozwjTGB2IbRHE1LKa/data/raw/ChIP-seq/TBL1XR1.bed?dl=1")
dest_locs <- paste0(offsite_dir, "/data/functional/ChIP-seq/", chip_seq_files, ".bed")
for (i in 1:length(chip_seq_files)) {
  download.file(url = urls[i], destfile = dest_locs[i])
}
