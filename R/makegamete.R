#' @title Make gamete
#'
#' @noRd
makegamete <- function(genotype) {
    gamete = genotype/2
    gamete[gamete==0.5] = sample(0:1, length(gamete[gamete==0.5]), replace=TRUE)
    return(gamete)
}
