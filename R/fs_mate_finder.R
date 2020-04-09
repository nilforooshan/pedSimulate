#' @title Find fullsib mates
#'
#' @noRd
fs_mate_finder <- function(mate, ped) {
  colnames(mate) = c("SIRE","DAM")
  colnames(ped) = c("ID","SIRE","DAM")
  mate = merge(mate, ped, by.x="SIRE", by.y="ID")
  colnames(mate) = c("SIRE","DAM","PGS","PGD")
  mate = merge(mate, ped, by.x="DAM", by.y="ID")
  colnames(mate) = c("DAM","SIRE","PGS","PGD","MGS","MGD")
  mate = mate[mate$PGS==mate$MGS & mate$PGD==mate$MGD, c("SIRE","DAM")]
  return(mate)
}
