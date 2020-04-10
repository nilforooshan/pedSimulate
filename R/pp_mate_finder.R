#' @title Find parent-progeny mates
#'
#' @noRd
pp_mate_finder <- function(mate, ped) {
  colnames(mate) = c("SIRE","DAM")
  colnames(ped) = c("ID","PGS","PGD")
  mate = merge(mate, ped[,c(1,3)], by.x="SIRE", by.y="ID")
  mate = mate[mate$PGD!=0,]
  colnames(ped) = c("ID","MGS","MGD")
  mate = merge(mate, ped[,1:2], by.x="DAM", by.y="ID")
  mate = mate[mate$MGS!=0,]
  mate = mate[mate$SIRE==mate$MGS | mate$DAM==mate$PGD, c("SIRE","DAM")]
  return(mate)
}
