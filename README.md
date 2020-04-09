# R package 'pedSimulate'

Simulate Pedigree, genetic merits and phenotypes with random mating followed by (non-)random selection differntly for males and females.

## Installation

```r
devtools::install_github('nilforooshan/pedSimulate')
```

## Description

An R package for simulating a pedigree with genetic merits and phenotypes, starting from a base population (generation 0).
The pedigree depth and design can be chosen by the values provided to the arguments of the simulation function.
Arguments are provided for the following: 

- Number of founder animals
- Proportion of females selected as dams
- Proportion of males selected as sires
- Mortality rate
- Litter size
- Number of generations to simulate
- Number of generation overlaps for sires
- Number of generation overlaps for dams
- (not )avoiding fullsib matings
- (not )avoiding parent-progeny matings
- Additive genetic variance in the base generation
- Environment (plus residual) variance
- (non-)random selection on females
- (non-)random selection on males

### For details, please read the PDF manual.

Thanks
