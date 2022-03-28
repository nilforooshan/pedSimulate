[![cran-version](https://www.r-pkg.org/badges/version/pedSimulate)](https://cran.r-project.org/package=pedSimulate)
[![cran-downloads-total](https://cranlogs.r-pkg.org/badges/grand-total/pedSimulate?color=green)](https://cran.r-project.org/package=pedSimulate)
[![cran-downloads-month](https://cranlogs.r-pkg.org/badges/last-month/pedSimulate?color=green)](https://cran.r-project.org/package=pedSimulate)
[![cran-downloads-week](https://cranlogs.r-pkg.org/badges/last-week/pedSimulate?color=green)](https://cran.r-project.org/package=pedSimulate)
[![github-release](https://img.shields.io/github/release/nilforooshan/pedSimulate.svg)](https://github.com/nilforooshan/pedSimulate)
[![github-license](https://img.shields.io/github/license/nilforooshan/pedSimulate.svg)](https://github.com/nilforooshan/pedSimulate/blob/master/LICENSE)
[![github-contributors](https://img.shields.io/github/contributors/nilforooshan/pedSimulate.svg)](https://github.com/nilforooshan/pedSimulate/graphs/contributors/)

# R package 'pedSimulate'

Simulate pedigree, genetic merits and phenotypes with random/assortative/disassortative matings followed by random/non-random selection of males and females with similar/different selection patterns in males and females.
For a given pedigree, genotypes can be simulated.

## Installation

Installing the package from CRAN:

```r
install.packages("pedSimulate")
```

Installing the package from GitHub:

```r
devtools::install_github('nilforooshan/pedSimulate')
```

## Description

An R package for simulating pedigree, genetic merits and phenotypes, starting from a base population (generation 0) or an existing pedigree.
The pedigree depth and design can be chosen by the values provided to the arguments of the simulation function.

Arguments are provided for the following: 

- Number of founder animals or an initial pedigree
- Additive genetic variance in the base generation
- Residual variance
- Litter size
- Number of generations to simulate after the founder generation
- Mortality rate
- Number of overlapping generations for sires
- Number of overlapping generations for dams
- Proportion of females selected as dams
- Proportion of males selected as sires
- Selection criterion for females
- Selection criterion for males
- Mating order for females
- Mating order for males

Genotypes can be simulated for a given pedigree. Marker effects to be chosen by the researcher.

Function `simulatePed` is for simulating a pedigree from a base population.

Function `appendPed` is for simulating new generations from an existing pedigree and appending to it.

Function `appendGen` is for simulating genotypes for an appended pedigree to an existing pedigree with genotypes.

Function `fs_mate_finder` is for finding fullsib matings in the pedigree.

Function `hs_mate_finder` is for finding halfsib matings in the pedigree.

Function `pp_mate_finder` is for finding parent-progeny matings in the pedigree.

Function `simulateGen` is for simulating genotypes for a given pedigree.

### For details, please read the PDF manual.

Thanks
