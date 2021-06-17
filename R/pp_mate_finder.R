#' @title Find parent-progeny mates
#'
#' @description Find parent-progeny matings in the pedigree
#'
#' @param ped : A pedigree \code{data.frame}. The first three columns (ID, SIRE, DAM) are used.
#'
#' @return pp_mates : A \code{data.frame} with two columns (SIRE, DAM) representing parent-progeny mates.
#'
#' @examples
#' ped = data.frame(ID=1:4, SIRE=c(0,0,1,1), DAM=c(0,0,2,3))
#' pp_mate_finder(ped)
#'
#' @export
pp_mate_finder <- function(ped) {
    ped = ped[,1:3]
    pp_mates = ped[,2:3]
    colnames(pp_mates) = c("SIRE","DAM")
    colnames(ped) = c("ID","PGS","PGD")
    pp_mates = pp_mates[pp_mates$SIRE!=0 & pp_mates$DAM!=0,]
    pp_mates = pp_mates[!duplicated(pp_mates),]
    pp_mates = merge(pp_mates, ped[,c(1,3)], by.x="SIRE", by.y="ID")
    colnames(ped) = c("ID","MGS","MGD")
    pp_mates = merge(pp_mates, ped[,1:2], by.x="DAM", by.y="ID")
    pp_mates = pp_mates[pp_mates$PGD!=0 | pp_mates$MGS!=0,]
    pp_mates = pp_mates[pp_mates$SIRE==pp_mates$MGS | pp_mates$DAM==pp_mates$PGD, c("SIRE","DAM")]
    return(pp_mates)
}
