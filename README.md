# SCEPTRE <img src="hex.png" align="right" width="200px"/>

SCEPTRE (analysis of **s**ingle **c**ell **p**er**t**urbation screens via conditional **re**sampling) is a method for single-cell CRISPR screen analysis. SCEPTRE infers gene-perturbation associations by modeling the stochastic assortment of CRISPR guide RNAs among cells instead of modeling gene expression, thereby remaining valid despite arbitrary misspecification of the gene expression model.

Accompanying paper:
> *Conditional resampling improves sensitivity and specificity of single cell CRISPR regulatory screens* <br />
> E. Katsevich, T. Barry, and K. Roeder (2020)<br />
> preprint available at [bioRxiv](https://doi.org/10.1101/2020.08.13.250092)


- [SCEPTRE methodology](#sceptre-methodology)
- [Repository contents and organization](#repository-contents-and-organization)
- [Dependencies](#dependencies)
- [Using this repository](#using-this-repository)
- [Code authors](#code-authors)

## SCEPTRE methodology

SCEPTRE proceeds one guide RNA (gRNA) and one gene at a time. It uses the z-value from a negative binomial regression to measure the effect of the gRNA on the gene. Instead of calibrating this z-value against a standard normal null distribution, it builds a null distribution for this statistic via conditional resampling. To this end, it first fits a logistic regression model for the occurrence of the gRNA in a cell, based on its covariates. For each cell, this yields a fitted probability that it contains the gRNA. Then, it generates a large number (default 500) of reshuffled datasets, where the expression and the covariates stay the same, while the gRNA assignment is redrawn independently for each cell based on its fitted probability. The negative binomial z-value is then recomputed for each of these datasets, which comprise a null distribution. A skew-t distribution is fit to these null histograms and the SCEPTRE p-value is defined by comparing the original z-value to this fitted skew-t null distribution.

<p align="center">
  <img src="sceptre_paper/manuscript/figures/Figure2/Figure2.png" width="600">
</p>

## Repository contents and organization

The repository contains two high-level directories: [sceptre_package](./sceptre_package) and [sceptre_paper](./sceptre_paper). The [sceptre_package](./sceptre_package) directory contains the `sceptre` R package and some shell scripts to help run `sceptre` at scale on computer clusters. The  [sceptre_paper](./sceptre_paper) directory contains code required to reproduce the analyses reported in Katsevich et al. 2020. Code in [sceptre_paper](./sceptre_paper) relies heavily on code in [sceptre_package](./sceptre_package); by contrast, code in [sceptre_package](./sceptre_package) does not depend at all on code in [sceptre_paper](./sceptre_paper).

* [sceptre_package](./sceptre_package): Contains the `sceptre` package and helper shell scripts.
  - [sceptre_package/sceptre](./sceptre_package/sceptre): `sceptre` package itself.
  - [sceptre_package/sceptre_at_scale](./sceptre_package/sceptre_at_scale): Helper shell scripts.
* [sceptre_paper](./sceptre_paper): Code for reproducing Katsevich et al. 2020.
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
- katsevich2020 0.1.0
- MASS 7.3.53
- monocle 2.18.0
- openxlsx 4.2.2
- R.matlab 3.6.2
- rhdf5 2.33.11
- sceptre 1.0.0
- Seurat 3.2.2
- sn 1.6.2
- tidyverse 1.3.0
- VGAM 1.1.3

Note that users looking to reproduce the analysis do not need to download these packages manually. Instead, the scripts to reproduce the analysis will (try to) download the packages automatically.

## Using this repository

Users can download and use the `sceptre` package independently of the manuscript analysis code. Note that the package currently is optimized to run the analyses reported in Katsevich et al. 2020. The authors are working on making the package more general and applicable to a wide range of datasets. If you would like to use SCEPTRE to analyze your data, please contact Eugene Katsevich at ekatsevi@wharton.upenn.edu.

#### Downloading and installing `sceptre`

Run the following code within R.

```
library(devtools)
install_github(repo="Timothy-Barry/SCEPTRE", subdir="sceptre_package/sceptre")
```

R may request that you update some packages on which `sceptre` relies. To be safe, please do so. If no dependencies need to be updated, then `sceptre` should download in a few seconds.

#### Reproducing the Katsevich et al. 2020 analysis

Git clone the SCEPTRE repository:

```
git clone https://github.com/Timothy-Barry/SCEPTRE.git
```

Next, navigate to the `sceptre_paper/utilities` directory, open the `run_everything.bash` script, and follow the instructions therein. The `run_everything.bash` script reproduces the entire analysis, from downloading the data and required packages to creating the figures. All analysis results (including intermediate files) are available on [Google Drive](https://drive.google.com/drive/folders/1ynZRMvGtFxfBiD0zAcuIYjNeS8Jj4AP9?usp=sharing). Parts of the analysis can therefore be reproduced by loading the intermediate results files instead of recomputing them.

#### Package demo and manual

A package demo is available [here](https://htmlpreview.github.io/?https://github.com/Timothy-Barry/SCEPTRE/blob/master/sceptre_package/sceptre/vignettes/sceptre-small-example.html). The code within the demo should take no longer than one minute to run. A package manual is available [here](https://github.com/Timothy-Barry/sceptre_paper/blob/master/sceptre_0.1.0.pdf).


## Code authors
- [Tim Barry](https://timothy-barry.github.io/)
- [Eugene Katsevich](https://ekatsevi.github.io/)
