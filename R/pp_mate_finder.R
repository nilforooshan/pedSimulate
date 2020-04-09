#' @title Find parent-progeny mates
#'
#' @noRd
pp_mate_finder <- function(mate, ped) {
  colnames(mate) = c("SIRE","DAM")
  colnames(ped) = c("ID","SIRE","DAM")
  mate = merge(mate, ped[,c(1,3)], by.x="SIRE", by.y="ID")
  colnames(mate) = c("SIRE","DAM","PGD")
  mate = merge(mate, ped[,1:2], by.x="DAM", by.y="ID")
  colnames(mate) = c("DAM","SIRE","PGD","MGS")
  mate = mate[mate$SIRE==mate$MGS | mate$DAM==mate$PGD, c("SIRE","DAM")]
  return(mate)
}
