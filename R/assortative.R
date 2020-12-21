#' @title Data simulation with assortative/disassortative matings
#'
#' @description Simulate Pedigree, genetic merits and phenotypes with assortative/disassortative mating
#' followed by non-random selection, differntly for males and females.
#'
#' @param F0size : Even number of founder animals. No mortality and selection in this generation, and matings are random in this generation.
#'
#' @param f.rate : Proportion of females selected as dams, default = 1.
#'
#' @param m.rate : Proportion of males (\code{<= f.rate}) selected as sires, default = 1.
#'
#' @param mort.rate : Mortality rate per generation, after the availability of phenotype (e.g., birth weight, weaning weight)
#' and before the age of maturity (i.e., before mating), default = 0. Maximum \code{mort.rate} = 0.5.
#'
#' @param littersize : Litter size, default = 1.
#'
#' @param ngen : Number of generations to simulate.
#'
#' @param overlap.s : Number of generation overlaps for sires, default = 0 for no generation overlap.
#'
#' @param overlap.d : Number of generation overlaps for dams, default = 0 for no generation overlap.
#'
#' @param Va0 : Additive genetic variance in the base generation (i.e., F0).
#'
#' @param Ve : Environment (plus residual) variance, set constant across generations.
#'
#' @param fsel : If \code{"P"} (default), selection on phenotypes or true breeding values if \code{Ve} = 0;
#' if \code{"PA"}, selection on true parent averages; redundant if \code{f.rate = 1}.
#' \code{"-P"} and \code{"-PA"} work in opposite direction of \code{"P"} and \code{"PA"}, respectively.
#'
#' @param msel : If \code{"P"} (default), selection on phenotypes or true breeding values if \code{Ve} = 0;
#' if \code{"PA"}, selection on true parent averages; redundant if \code{m.rate = 1}.
#' \code{"-P"} and \code{"-PA"} work in opposite direction of \code{"P"} and \code{"PA"}, respectively.
#'
#' @param negative : Assortative mating if \code{TRUE} (default) and disassortative mating if \code{false}.
#' Males are sorted based on \code{msel} and females are sorted based on \code{fsel}.
#'
#' @return ped : The output pedigree \code{data.frame}. Further information provided in \strong{Details}.
#'
#' @details
#' The output pedigree \code{data.frame} (\code{ped}) has 9 columns: ID, SIRE, DAM, SEX,
#' GEN (generation number starting with 0 for the base generation), PA (parent average),
#' MS (Mendelian Sampling), E (environment and residuals), and P (phenotype).
#'
#' @examples
#' ped = assortative(
#'     F0size = 100,
#'     f.rate = 0.8,
#'     m.rate = 0.5,
#'     mort.rate = 0.05,
#'     littersize = 2,
#'     ngen = 4,
#'     overlap.s = 1,
#'     overlap.d = 0,
#'     negative = FALSE,
#'     Va0 = 9,
#'     Ve = 36,
#'     fsel = "P",
#'     msel = "PA"
#' )
#'
#' @export
assortative <- function(F0size, f.rate=1, m.rate=1, mort.rate=0, littersize=1, ngen, overlap.s=0, overlap.d=0, Va0, Ve, fsel="P", msel="P", negative=FALSE) {
    # Check inputs
    F0size = round(F0size)
    littersize = round(littersize)
    ngen = round(ngen)
    overlap.s = round(overlap.s)
    overlap.d = round(overlap.d)
    if(F0size < 2) stop("ERROR: F0size < 2")
    if((F0size %% 2)!=0) stop("ERROR: F0size should be an even number.")
    if(f.rate > 1) stop("ERROR: f.rate > 1")
    if(f.rate <= 0) stop("ERROR: f.rate <= 0")
    if(m.rate > 1) stop("ERROR: m.rate > 1")
    if(m.rate <= 0) stop("ERROR: m.rate <= 0")
    if(mort.rate > 0.5) stop("ERROR: mort.rate > 0.5")
    if(mort.rate < 0) stop("ERROR: mort.rate < 0")
    if(f.rate < m.rate) stop("ERROR: f.rate < m.rate")
    if(littersize < 1) stop("ERROR: littersize < 1")
    if(ngen < 1) stop("ERROR: ngen < 1")
    if(overlap.s < 0) stop("ERROR: overlap.s < 0")
    if(overlap.d < 0) stop("ERROR: overlap.d < 0")
    if(overlap.s >= ngen) stop("ERROR: overlap.s >= ngen")
    if(overlap.d >= ngen) stop("ERROR: overlap.d >= ngen")
    if(Va0 < 0) stop("ERROR: Va0 < 0")
    if(Ve < 0) stop("ERROR: Ve < 0")
    if(!fsel %in% c("P","PA","-P","-PA")) stop('ERROR: fsel should be "P", "PA", "-P" or "-PA".')
    if(!msel %in% c("P","PA","-P","-PA")) stop('ERROR: msel should be "P", "PA", "-P" or "-PA".')
    if(!negative %in% c(TRUE,FALSE)) stop("ERROR: negative should be TRUE or FALSE.")
    # F0
    ordcol = c("ID","SIRE","DAM","SEX","GEN","SBV","DBV","MS","E","P","DEAD")
    ped = data.frame(ID=1:F0size,
                     SIRE=0,
                     DAM=0,
                     SEX=rep(c("m","f"), each=F0size/2),
                     GEN=0,
                     SBV=0,
                     DBV=0,
                     MS=rnorm(F0size, 0, sqrt(Va0)),
                     E=rnorm(F0size, 0, sqrt(Ve)))
    ped$P = ped$MS + ped$E
    ped$DEAD = FALSE
    message("Started with a base generation of ", F0size, " individuals.")
    # F1
    tmp = data.frame(ID=(1:(F0size*littersize/2))+F0size,
                     SIRE=rep(ped[ped$SEX=="m",]$ID, littersize),
                     DAM=rep(sample(ped[ped$SEX=="f",]$ID), each=littersize))
    tmp = merge(tmp, ped[,c("ID","SBV","DBV","MS")], by.x="SIRE", by.y="ID")
    tmp = merge(tmp, ped[,c("ID","SBV","DBV","MS")], by.x="DAM",  by.y="ID")
    tmp$SBV = (tmp$SBV.x + tmp$DBV.x)/2 + tmp$MS.x
    tmp$DBV = (tmp$SBV.y + tmp$DBV.y)/2 + tmp$MS.y
    tmp$SEX= sample(c("m","f"), F0size*littersize/2, replace=TRUE)
    tmp$GEN=1
    tmp$MS=rnorm(nrow(tmp), 0, sqrt(Va0/2))
    tmp$E=rnorm(nrow(tmp), 0, sqrt(Ve))
    tmp$P = (tmp$SBV + tmp$DBV)/2 + tmp$MS + tmp$E
    tmp$DEAD = FALSE
    message("Simulated generation 1 with ", nrow(tmp), " individuals.")
    ped = rbind(ped, tmp[,ordcol])
    # Fi (i>1)
    if(ngen > 1) {
        for(i in 2:ngen)
        {
            # Mortality before maturity
            if(mort.rate > 0) {
                ndead = round(mort.rate*nrow(ped[ped$GEN==(i-1),]))
                ped[sample(which(ped$GEN==(i-1)), ndead), "DEAD"] = TRUE
            }
            sires = ped[ped$SEX=="m" & ped$GEN %in% (-overlap.s:0)+i-1 & !ped$DEAD, c("ID","SBV","DBV","P")]
            dams  = ped[ped$SEX=="f" & ped$GEN %in% (-overlap.d:0)+i-1 & !ped$DEAD, c("ID","SBV","DBV","P")]
            # Selection on males
            nm = round(m.rate*nrow(sires))
            if(nm==0) stop("ERROR: No male left. Generation ", i)
            if(msel=="P") {
                sires = sires[order(-sires$P),]$ID[1:nm]
            } else if(msel=="PA") {
                sires = sires[order(-sires$SBV -sires$DBV),]$ID[1:nm]
            } else if(msel=="-P") {
                sires = sires[order(sires$P),]$ID[1:nm]
            } else if(msel=="-PA") {
                sires = sires[order(sires$SBV + sires$DBV),]$ID[1:nm]
            }
            # Selection on females
            if(negative) {
                pn = 1
            } else {
                pn = -1
            }
            nf = round(f.rate*nrow(dams))
            if(fsel=="P") {
                dams = dams[order(pn*dams$P),]$ID[1:nf]
            } else if(fsel=="PA") {
                dams = dams[order(pn*dams$SBV + pn*dams$DBV),]$ID[1:nf]
            } else if(fsel=="-P") {
                dams = dams[order(-pn*dams$P),]$ID[1:nf]
            } else if(fsel=="-PA") {
                dams = dams[order(-pn*dams$SBV -pn*dams$DBV),]$ID[1:nf]
            }
            # Set mates
            SIRES = c(rep(sires, each=floor(nf/nm)), sires[1:(nf-length(sires)*floor(nf/nm))])
            SIRES = SIRES[order(match(SIRES, sires))]
            tmp = data.frame(SIRE=SIRES, DAM=dams)
            # Append the next generation
            tmp = data.frame(SIRE=rep(tmp$SIRE, littersize), DAM=rep(tmp$DAM, littersize))
            tmp$ID = (1:nrow(tmp))+nrow(ped)
            tmp = merge(tmp, ped[,c("ID","SBV","DBV","MS")], by.x="SIRE", by.y="ID")
            tmp = merge(tmp, ped[,c("ID","SBV","DBV","MS")], by.x="DAM",  by.y="ID")
            tmp$SBV = (tmp$SBV.x + tmp$DBV.x)/2 + tmp$MS.x
            tmp$DBV = (tmp$SBV.y + tmp$DBV.y)/2 + tmp$MS.y
            tmp$SEX = sample(c("m","f"), nrow(tmp), replace=TRUE)
            tmp$GEN = i
            tmp$MS=rnorm(nrow(tmp), 0, sqrt(Va0/2))
            tmp$E=rnorm(nrow(tmp), 0, sqrt(Ve))
            tmp$P = (tmp$SBV + tmp$DBV)/2 + tmp$MS + tmp$E
            tmp$DEAD = FALSE
            message("Simulated generation ", i, " with ", nrow(tmp), " individuals.")
            ped = rbind(ped, tmp[,ordcol])
        }
    }
    # Replace SBV & DBV with PA
    ped$PA = (ped$SBV + ped$DBV)/2
    ped = ped[,c("ID","SIRE","DAM","SEX","GEN","PA","MS","E","P")]
    # Re-base P
    ped$P = ped$P - 2*min(ped$P)
    ped = ped[order(ped$ID),]
    rownames(ped) = 1:nrow(ped)
    return(ped)
}
