#' @title Find fullsib mates
#'
#' @description Find fullsib matings in the pedigree
#'
#' @param ped : A pedigree \code{data.frame}. The first three columns (ID, SIRE, DAM) are used.
#'
#' @return fs_mates : A \code{data.frame} with two columns (SIRE, DAM) representing fullsib mates.
#'
#' @examples
#' ped = data.frame(ID=1:7, SIRE=c(0,0,1,0,3,3,5), DAM=c(0,0,0,2,4,4,6))
#' fs_mate_finder(ped)
#'
#' @export
fs_mate_finder <- function(ped) {
    ped = ped[,1:3]
    fs_mates = ped[,2:3]
    colnames(fs_mates) = c("SIRE","DAM")
    colnames(ped) = c("ID","PGS","PGD")
    fs_mates = fs_mates[fs_mates$SIRE!=0 & fs_mates$DAM!=0,]
    fs_mates = fs_mates[!duplicated(fs_mates),]
    fs_mates = merge(fs_mates, ped, by.x="SIRE", by.y="ID")
    fs_mates = fs_mates[fs_mates$PGS!=0 & fs_mates$PGD!=0,]
    colnames(ped) = c("ID","MGS","MGD")
    fs_mates = merge(fs_mates, ped, by.x="DAM", by.y="ID")
    fs_mates = fs_mates[fs_mates$MGS!=0 & fs_mates$MGD!=0,]
    fs_mates = fs_mates[fs_mates$PGS==fs_mates$MGS & fs_mates$PGD==fs_mates$MGD, c("SIRE","DAM")]
    return(fs_mates)
}
