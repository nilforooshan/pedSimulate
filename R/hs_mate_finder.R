#' @title Find halfsib mates
#'
#' @description Find halfsib matings in the pedigree
#'
#' @param ped : A pedigree \code{data.frame}. The first three columns (ID, SIRE, DAM) are used.
#'
#' @return hs_mates : A \code{data.frame} with two columns (SIRE, DAM) representing halfsib mates.
#'
#' @examples
#' ped = data.frame(ID=1:7, SIRE=c(0,0,1,1,0,3,5), DAM=c(0,0,2,2,2,4,4))
#' hs_mate_finder(ped)
#'
#' @export
hs_mate_finder <- function(ped) {
    ped = ped[,1:3]
    hs_mates = ped[,2:3]
    colnames(hs_mates) = c("SIRE","DAM")
    colnames(ped) = c("ID","PGS","PGD")
    hs_mates = hs_mates[hs_mates$SIRE!=0 & hs_mates$DAM!=0,]
    hs_mates = hs_mates[!duplicated(hs_mates),]
    hs_mates = merge(hs_mates, ped, by.x="SIRE", by.y="ID")
    hs_mates = hs_mates[hs_mates$PGS!=0 | hs_mates$PGD!=0,]
    colnames(ped) = c("ID","MGS","MGD")
    hs_mates = merge(hs_mates, ped, by.x="DAM", by.y="ID")
    hs_mates = hs_mates[hs_mates$MGS!=0 | hs_mates$MGD!=0,]
    hs_mates = hs_mates[(hs_mates$MGS!=0 & hs_mates$MGS==hs_mates$PGS) | (hs_mates$MGD!=0 & hs_mates$MGD==hs_mates$PGD),]
    # Discard fullsibs
    hs_mates = hs_mates[hs_mates$MGS==0 | hs_mates$MGD==0 | hs_mates$MGS!=hs_mates$PGS | hs_mates$MGD!=hs_mates$PGD, c("SIRE","DAM")]
    return(hs_mates)
}
