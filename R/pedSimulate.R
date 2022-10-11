#' @title Pedigree, genetic merit, phenotype, and genotype simulation
#'
#' @docType package
#'
#' @name pedSimulte-package
#'
#' @author Mohammad Ali Nilforooshan \email{m.a.nilforooshan@gmail.com}
#'
#' @description
#' An R package for simulating a pedigree with genetic merits and phenotypes, starting from a base population (generation 0) or an existing pedigree.
#' The pedigree depth and design can be chosen by the values provided to the arguments of the simulation function.
#' Genotypes can be simulated for a given pedigree, or an appended pedigree to an existing pedigree with genotypes.
#' Marker effects to be chosen by the researcher.
#'
#' @details
#' Starting from a base population with a user-defined size and equal number of males and females,
#' next generations are simulated for the user-defined litter size and number of generations.
#' No selection (natural or artificial) and non-random mating is applied to this population.
#' Alternatively, the simulation can be started from an existing pedigree.
#' Natural (mortality) and artificial selection are applied to the next generations.
#' Different generation overlap, selection intensities and selection patterns can be applied to males and females.
#' Selected males and females are ordered similarly/differently to simulate various random, assortative or disassortative mating scenarios.
#' Performance and genetic merit of individuals are simulated using the basic rules of quantitative genetics.
#' The performance (P) of an individual is influenced by genetic (A) and environmental (E) effects.
#' Thus, P = A + E, and Var(P) = Var(A) + Var(E).
#' The additive genetic merit (A) of an individual is the average of its parents' additive genetic merits
#' (PA = (A\ifelse{html}{\out{<sub>sire</sub>}}{\eqn{_{sire}}} + A\ifelse{html}{\out{<sub>dam</sub>}}{\eqn{_{dam}}})/2)
#' plus the Mendelian Sampling term due to the sampling of alleles passed from the parent to the offspring.
#' The Mendelian Sampling variance is half of Var(A) in the base population.
#' Because there is no provided information for environmental effects, the environment effect is
#' assigned to individuals from a normal distribution of random numbers (E ~ N(0, \strong{I}Var(E))).
#' The package also provides functions to identify halfsib, fullsib and parent-progeny matings in the pedigree.
#' For a given pedigree, genotypes can be simulated.
#' Marker effects can be chosen by the researcher and added to the simulated phenotypes.
#' In that case, genetic effects and variance used to simulate phenotypes change to residual polygenic effects and variance (genetic effects and variance not explained by the markers).
#'
#' @references
#' Mrode, R. A. 2005. Linear Models for the Prediction of Animal Breeding Values, 2nd ed. Cambridge, MA: CABI Publishing. <ISBN:9780851989969, 0851989969>
#'
#' Nilforooshan, M. A. 2022. pedSimulate – An R package for simulating pedigree, genetic merit, phenotype, and genotype data. R. Bras. Zootec., 51:e20210131. <doi:10.37496/rbz5120210131>
NULL
