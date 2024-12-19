2024-12-18

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication repository for “Forecasting macroeconomic data with Bayesian VARs: Sparse or dense? It depends!”

Luis Gruber (<luis.gruber@aau.at>) and Gregor Kastner.

Preprint available at
[Arxiv](https://doi.org/10.48550/arXiv.2206.04902). Working paper is
conditionally accepted at [International Journal of
Forecasting](https://forecasters.org/ijf/).

## Overview

This repository provides the replication material for the paper
“Forecasting macroeconomic data with Bayesian VARs: Sparse or dense? It
depends!”. Each figure and table is generated separately by its
corresponding script file `figure_[x].R` or `table_[x].R`, respectively.

The main contents of the repository are the following:

- `tables/`: Folder containing the generated tables as .txt files.
- `figures/`: Folder containing the generated figures as .pdf or .png
  files.
- `data/`: Folder containing the data used in the empirical application.
- `results/`: Folder containing intermediary results which were
  generated using a computer cluster.
- `computer_cluster_scripts/`: Folder containing the computer cluster
  scripts.
- `figure_[x].R`: R scripts to generate the respective figures.
- `table_[x].R`: R scripts to generate the respective tables.
- `bayesianVARs_0.1.3.tar.gz`: Source files of R package
  [bayesianVARs](https://luisgruber.github.io/bayesianVARs/) (0.1.3).
- `bayesianVARs_0.1.3.zip`: `.zip` file containing the binary builds of
  the R package
  [bayesianVARs](https://luisgruber.github.io/bayesianVARs/) (0.1.3).

## Instruction & computational requirements

All file paths are relative to the root of the replication repository.
Please set your working directory accordingly.

All the estimation and analysis is done in R. The main package used for
estimation and forecasting is `bayesianVARs` (0.1.3) with major
dependency `stochvol` (3.2.0), where the values in parenthesis indicate
the package versions we used. In addition, the following add-on packages
are loaded and attached in the replication scripts (in alphabetical
order): `coda` (0.19-4), `colorspace` (2.1-1), `dplyr` (1.1.4),
`ggplot2` (3.5.1), `gsl` (2.1-8), `knitr` (1.49), `lubridate` (1.9.3),
`mvtnorm` (1.3-2), `purrr` (1.0.2), `readr` (2.1.5), `reticulate`
(1.40.0), `RhpcBLASctl` (0.23-42), `scatterplot3d` (0.3-44), `stringr`
(1.5.1), `tibble` (3.2.1), `tidyr` (1.3.1) and `xts` (0.13.0).

The files `figure_[x].R` and `table_[x].R` can be run individually in
any order. On a standard computer the execution of each those files
should take less than one minute, except `figure_3.R`, which should take
less than five minutes.

We provide the source code of `bayesianVARs` in the file
[bayesianVARs_0.1.3.tar.gz](bayesianVARs_0.1.3.tar.gz). A great deal of
the package is written in C++. Therefore, suitable compilers and related
tools (e.g. Rtools on Windows) need to be installed in order to install
the package from source. Additionally, we provide the file
[bayesianVARs_0.1.3.zip](bayesianVARs_0.1.3.zip) containing the binary
builds. E.g., on Windows you can install the package in R with
`install.packages("bayesianVARs_0.1.3.zip")`.

The Python library `mpmath` is required in order to fully reproduce
[figures/figure_1.png](./figures/figure_1.png). Via R this can be done
through the package `reticulate` (see the instructions in
[figure_1.R](figure_1.R)).

The computations of the main results, more concretely the simulation
study in Section 4.2 and the empirical application in Section 5, were
carried out on a computer cluster with 3328 CPUs managed with the [Slurm
Workload Manager](https://slurm.schedmd.com/). The corresponding R
scripts are located in the folder
[computer_cluster_scripts/](./computer_cluster_scripts/). Those scripts
should not be executed on standard computers. The script
[computer_cluster_scripts/simulation_study_clusterscript.R](./computer_cluster_scripts/simulation_study_clusterscript.R)
estimates 3240 models. The script
[computer_cluster_scripts/empirical_application_clusterscript.R](./computer_cluster_scripts/empirical_application_clusterscript.R)
estimates and evaluates forecasts of 284970 models. Even on a computer
cluster it can take several days or even weeks to run those scripts,
depending on the available resources.

## Data

The data of the empirical application consists of 21 selected quarterly
time-series obtained from the FRED-QD data base (McCracken and Ng, 2021,
<https://doi.org/10.20955/r.103.1-44>). Release date 2021-07
(<https://files.stlouisfed.org/files/htdocs/fred-md/quarterly/2021-07.csv>).
The data is transformed to be interpreted as growth-rates (first
log-differences with the exception of interest rates, which are already
growth rates). The prepared data is located at
[data/data_growth.RData](./data/data_growth.RData).

In Figures 6 through 8 the NBER recession are highlighted among other
things. We obtained the relevant data from the R package
[tis](https://cran.r-project.org/package=tis) and provide it as `.rds`
file at [data/nber_recessions.rds](./data/nber_recessions.rds).
