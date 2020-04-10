#' @title Find fullsib mates
#'
#' @noRd
fs_mate_finder <- function(mate, ped) {
  colnames(mate) = c("SIRE","DAM")
  colnames(ped) = c("ID","PGS","PGD")
  mate = merge(mate, ped, by.x="SIRE", by.y="ID")
  mate = mate[mate$PGS!=0 & mate$PGD!=0,]
  colnames(ped) = c("ID","MGS","MGD")
  mate = merge(mate, ped, by.x="DAM", by.y="ID")
  mate = mate[mate$MGS!=0 & mate$MGD!=0,]
  mate = mate[mate$PGS==mate$MGS & mate$PGD==mate$MGD, c("SIRE","DAM")]
  return(mate)
}
