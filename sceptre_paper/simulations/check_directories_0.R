# Check directory structure
args <- commandArgs(trailingOnly = TRUE)
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]

sub_dirs <- c("data/simulations", "results/simulations", "logs/simulations")
dirs_to_create <- paste0(offsite_dir, "/", sub_dirs)
for (directory in dirs_to_create) {
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
}
