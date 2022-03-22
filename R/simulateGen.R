#' @title Simulate genotypes
#'
#' @description Simulate genotypes for a given pedigree, allele frequency and mutation rate at each marker locus.
#'
#' @param ped : Pedigree \code{data.frame} with columns for animal, sire, and dam identification.
#'
#' @param AF : Vector of allele frequencies at different loci for the genotypes to be simulated.
#'
#' @param mut.rate : Vector of mutation rates at different loci for the genotypes to be simulated, default = 0 for no mutation.
#'
#' @param seed : A numeric variable input to the random number generator for reproducible simulations,
#' default = `NA` for non-reproducible simulations.
#'
#' @return M : The simulated genotype \code{data.frame} with rows corresponding to animals (in the same order as in the pedigree) and columns corresponding to markers.
#'
#' @details
#' Only diploid and bi-allelic situations are covered.
#' No linkage disequilibrium is simulated.
#'
#' @examples
#' nSNP = 100
#' AF = runif(nSNP, 0.01, 0.99)
#' mut.rate = runif(nSNP, 0, 10^-5)
#' ped = data.frame(ID=1:5, SIRE=c(0,0,1,0,3), DAM=c(0,0,2,2,4))
#' gen = simulateGen(ped, AF, mut.rate, seed=684)
#'
#' @export
simulateGen <- function(ped, AF, mut.rate=0, seed=NA) {
    colnames(ped) = c("ID","SIRE","DAM")
    stopifnot(all(is.numeric(AF)))
    stopifnot(all(!is.na(AF)))
    stopifnot(min(AF)>=0.01)
    stopifnot(max(AF)<=0.99)
    if(!identical(mut.rate, 0)) {
        stopifnot(length(AF)==length(mut.rate))
        stopifnot(min(mut.rate)>=0)
        if(length(mut.rate[mut.rate > 10^-6])) {
            warning("Found ", length(mut.rate[mut.rate > 10^-6]), " markers with mutation rate > 10^-6")
            warning("Maximum mutation rate = ", max(mut.rate))
        }
    } else {
        message("No mutation was simulated.")
    }
    stopifnot(length(seed)==1)
    stopifnot(is.na(seed) | is.numeric(seed))
    if(!is.na(seed)) set.seed(seed)
    tmp = c()
    SNPs = 1:length(AF)
    # Imputing the 1st genotype
    for(j in SNPs) tmp = c(tmp, sample(0:2, 1, prob=c((1-AF[j])^2, 2*(1-AF[j])*AF[j], AF[j]^2)))
    M = matrix(tmp, nrow=1)
    # Imputing the next genotypes
    for(i in 2:nrow(ped))
    {
        s = ped$SIRE[i]
        d = ped$DAM[i]
        tmp = c()
        if(s==0 & d==0) {
            for(j in SNPs) tmp = c(tmp, sample(0:2, 1, prob=c((1-AF[j])^2, 2*(1-AF[j])*AF[j], AF[j]^2)))
        } else if(s>0 & d==0) {
            for(j in SNPs) tmp = c(tmp, sample(0:1, 1, prob=c(1-AF[j], AF[j])))
            tmp = tmp+makegamete(M[s,], seed=seed)
        } else if(s==0 & d>0) {
            for(j in SNPs) tmp = c(tmp, sample(0:1, 1, prob=c(1-AF[j], AF[j])))
            tmp = tmp+makegamete(M[d,], seed=seed)
        } else {
            tmp = makegamete(M[s,], seed=seed)+makegamete(M[d,], seed=seed)
        }
        if(!identical(mut.rate, 0)) tmp = mutate(tmp, mut.rate, seed=seed)
        M = rbind(M, tmp)
        rownames(M) = NULL
    }
    return(M)
}
