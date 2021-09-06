#' @title Make mutation
#'
#' @noRd
mutate <- function(genotype, mut.rate) {
    stopifnot(length(genotype)==length(mut.rate))
    for(i in which(mut.rate > 0))
    {
        if(genotype[i]==0) {
            genotype[i] = sample(0:2, 1, prob=c(1-mut.rate[i]-mut.rate[i]^2, mut.rate[i], mut.rate[i]^2))
        } else if(genotype[i]==1) {
            genotype[i] = sample(0:2, 1, prob=c(mut.rate[i], 1-2*mut.rate[i], mut.rate[i]))
        } else {
            genotype[i] = sample(0:2, 1, prob=c(mut.rate[i]^2, mut.rate[i], 1-mut.rate[i]-mut.rate[i]^2))
        }
    }
    return(genotype)
}
