#' @title Pedigree, genetic merit and phenotype simulation
#'
#' @docType package
#'
#' @name pedSimulte-package
#'
#' @author Mohammad Ali Nilforooshan \email{m.a.nilforooshan@gmail.com}
#'
#' @description
#' An R package for simulating a pedigree with genetic merits and phenotypes, starting from a base population (generation 0).
#' The pedigree depth and design can be chosen by the values provided to the arguments of the simulation function.
#'
#' @details
#' Starting from a base population with equal number of males and females, next generations are simulated for a user-defined number of generations.
#' The size of the base population is defined by the user, and there is no natural or artificial selection applied to this population.
#' Natural (mortality) and artificial selection are applied to the next generations.
#' Different selection patterns can be applied to males and females, including
#' the proportion of selected individuals from selection candidates (after applying mortality) as parents of the next generation,
#' random or merit-based selection,
#' and generation overlap.
#' Selected individuals are mated randomly.
#' Further choices, such as litter size, and avoiding fullsib matings and parent-progeny matings are provided.
#' Performance and genetic merit of individuals are simulated using the basic rules of quantitative genetics.
#' The performance (P) of an individual is influenced by genetic (A) and environmental (E) effects.
#' Thus, P = A + E, and Var(P) = Var(A) + Var(E).
#' The additive genetic merit (A) of an individual is the average of its parents' additive genetic merits (PA = (A\ifelse{html}{\out{<sub>sire</sub>}}{\eqn{_{sire}}} + A\ifelse{html}{\out{<sub>dam</sub>}}{\eqn{_{dam}}})/2)
#' plus the Mendelian Sampling term (MS) due to the sampling of alleles passed from the parent to the offspring.
#' The variance of MS is equal to half of Var(A) in the base population.
#' Because there is no provided information for environmental effects, the environment effect is assigned to individuals from a normal distribution of random numbers (E ~ N(0, \strong{I}Var(E))).
NULL
