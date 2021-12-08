#' @title Simulate genotypes for an appended pedigree
#'
#' @description Simulate genotypes for an appended pedigree to an existing pedigree with genotypes.
#'
#' @param ped : Pedigree \code{data.frame} with columns for animal, sire, and dam identification.
#'
#' @param M : Genotype \code{data.frame} with rows corresponding to the initial rows of the pedigree and columns corresponding to markers.
#'
#' @param AF : Vector of allele frequencies at different loci for the genotypes to be simulated. If no value is provided, it will be estimated from \code{M}.
#'
#' @param mut.rate : Vector of mutation rates at different loci for the genotypes to be simulated, default = 0 for no mutation.
#'
#' @return M2 : New simulated genotypes appended to \code{M}.
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
#' gen = simulateGen(ped, AF, mut.rate)
#' ped = rbind(ped, data.frame(ID=6:8, SIRE=c(3,6,6), DAM=c(0,4,5)))
#' gen = appendGen(ped, M=gen, AF)
#'
#' @export
appendGen <- function(ped, M, AF=c(), mut.rate=0) {
    if(!identical(as.integer(ped[,1]), 1:nrow(ped))) stop("!identical(as.integer(ped[,1]), 1:nrow(ped)); Please consider renumberring the pedigree: ped[,1:3] <- ggroups::renum(ped[,1:3])$newped")
    stopifnot(nrow(M) < nrow(ped))
    if(nrow(M)==0) stop("The genotype matrix is empty. Please use function simulateGen.")
    colnames(ped) = c("ID","SIRE","DAM")
    if(length(AF)==0) {
        AF = colMeans(M)
    } else {
        stopifnot(min(AF)>=0.01)
        stopifnot(max(AF)<=0.99)
    }
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
    tmp = c()
    SNPs = 1:length(AF)
    # Appending genotypes
    for(i in (nrow(M)+1):nrow(ped))
    {
        s = ped$SIRE[i]
        d = ped$DAM[i]
        tmp = c()
        if(s==0 & d==0) {
            for(j in SNPs) tmp = c(tmp, sample(0:2, 1, prob=c((1-AF[j])^2, 2*(1-AF[j])*AF[j], AF[j]^2)))
        } else if(s>0 & d==0) {
            for(j in SNPs) tmp = c(tmp, sample(0:1, 1, prob=c(1-AF[j], AF[j])))
            tmp = tmp+makegamete(M[s,])
        } else if(s==0 & d>0) {
            for(j in SNPs) tmp = c(tmp, sample(0:1, 1, prob=c(1-AF[j], AF[j])))
            tmp = tmp+makegamete(M[d,])
        } else {
            tmp = makegamete(M[s,])+makegamete(M[d,])
        }
        if(!identical(mut.rate, 0)) tmp = mutate(tmp, mut.rate)
        M = rbind(M, tmp)
        rownames(M) = NULL
    }
    return(M)
}
