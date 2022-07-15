# SCEPTRE

SCEPTRE is a method for single-cell CRISPR screen analysis. This repository contains code to reproduce all analyses reported in the following paper, which introduces the SCEPTRE method:

> *Conditional resampling improves calibration and sensitivity in single-cell CRISPR screen analysis* <br />
> T. Barry, X. Wang, J. Morris, K. Roeder, and E. Katsevich (2021)<br />
> preprint available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.08.13.250092v6)

Note that there are two `sceptre` R packages: an in-house package intended for manuscript reproduction purposes only, and a separate package intended for new data analysis purposes. If you would like to use SCEPTRE to analyze your own data, please visit the [Github repository](https://github.com/Katsevich-Lab/sceptre) or [website](https://katsevich-lab.github.io/sceptre/) of the latter package.

## Repository contents and organization

This repository contains two high-level directories: [sceptre_package](./sceptre_package) and [sceptre_paper](./sceptre_paper). The [sceptre_package](./sceptre_package) directory contains the (in-house) `sceptre` R package and some shell and R scripts to help run the method at scale. The [sceptre_paper](./sceptre_paper) directory contains code required to reproduce all analyses reported in Barry et al. 2020. Code in [sceptre_paper](./sceptre_paper) relies heavily on code in [sceptre_package](./sceptre_package); by contrast, code in [sceptre_package](./sceptre_package) does not depend at all on code in [sceptre_paper](./sceptre_paper).

* [sceptre_package](./sceptre_package): Contains the (in-house) `sceptre` package and helper scripts.
  - [sceptre_package/sceptre](./sceptre_package/sceptre): (in-house) `sceptre` package itself.
  - [sceptre_package/sceptre_at_scale](./sceptre_package/sceptre_at_scale): Helper shell and R scripts.
* [sceptre_paper](./sceptre_paper): Code for reproducing Barry et al. 2020.
  - [sceptre_paper/analysis_drivers](./sceptre_paper/analysis_drivers): R scripts for reproducing data analysis.
  - [sceptre_paper/katsevich2020](./sceptre_paper/katsevich2020): Additional R package containing functions called by scripts in [sceptre_paper](./sceptre_paper) directory.
  - [sceptre_paper/manuscript](./sceptre_paper/manuscript): .tex manuscript.
  - [sceptre_paper/nb_regression_at_scale](./sceptre_paper/nb_regression_at_scale): Code for running large-scale negative binomial regression.
  - [sceptre_paper/plotting](./sceptre_paper/plotting): R scripts to create manuscript figures.
  - [sceptre_paper/simulations](./sceptre_paper/simulations): R scripts for running simulations.
  - [sceptre_paper/utilities](./sceptre_paper/utilities): Shell and R scripts for reproducing the entire analysis.

## Dependencies

The code was executed across two machines: a computer cluster running R version 3.6.1 and Linux kernel 4.4.180-102, and a Macbook running R version 4.0.2 and Darwin kernel version 20.1.0. In addition, the following R packages were used:

- bigstatsr 1.2.3
- cowplot 1.1.0
- fst 0.9.4
- furrr 0.2.0
- ggpubr 0.4.0
- ggrepel 0.9.1
- katsevich2020 0.1.0
- MASS 7.3.53
- mgcv 1.8-33
- monocle 2.18.0
- openxlsx 4.2.2
- R.matlab 3.6.2
- rhdf5 2.33.11
- scales 1.1.1
- sceptre 1.0.0
- Seurat 3.2.2
- sn 1.6.2
- tidyverse 1.3.0
- VGAM 1.1.3

The scripts to reproduce the analysis will (try to) download the packages automatically.

#### Reproducing the Barry et al. 2020 analysis

Git clone the SCEPTRE repository:

```
git clone https://github.com/Katsevich-Lab/sceptre-manuscript
```

Next, navigate to the `sceptre_paper/utilities` directory, open the `run_everything.bash` script, and follow the instructions therein. The `run_everything.bash` script reproduces the entire analysis, from downloading the data to creating the figures. All analysis results (including intermediate files) are available on [Box](https://upenn.box.com/v/sceptre-files-v8). Parts of the analysis therefore can be reproduced by loading the intermediate results files instead of recomputing them.

## Contributions

### Code authors
- [Tim Barry](https://timothy-barry.github.io/)
- [Eugene Katsevich](https://ekatsevi.github.io/)
- [Xuran Wang](https://xuranw.github.io/personalwebsite/)

### Code contributors
- [John A. Morris](https://morrisgenetics.com/publications/)
