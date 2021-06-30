#' @title Simulate pedigree, genetic merits and phenotypes
#'
#' @description Simulate pedigree, genetic merits and phenotypes with random/assortative/disassortative matings
#' followed by random/non-random selection of males and females with similar/different patterns in males and females.
#'
#' @param F0size : Even number of founder animals. No mortality, selection and non-random mating in this generation.
#'
#' @param Va0 : Additive genetic variance in the base generation (i.e., F0).
#'
#' @param Ve : Residual variance, constant across generations.
#'
#' @param littersize : Litter size, default = 1.
#'
#' @param ngen : Number of generations to simulate.
#'
#' @param mort.rate : Mortality rate per generation, after the availability of phenotype (e.g., birth weight, weaning weight)
#' and before the age of maturity (i.e., before mating), default = 0. Maximum \code{mort.rate} = 0.5.
#'
#' @param overlap.s : Number of overlapping generations for sires, default = 0 for no generation overlap.
#'
#' @param overlap.d : Number of overlapping generations for dams, default = 0 for no generation overlap.
#'
#' @param f.rate : Proportion of females selected as dams, default = 1.
#'
#' @param m.rate : Proportion of males (\code{<= f.rate}) selected as sires, default = 1.
#'
#' @param fsel : If \code{"R"} (default), random selection on females;
#' if \code{"P"}, selection on phenotypes or true breeding values if \code{Ve} = 0;
#' if \code{"PA"}, selection on true parent averages.
#' \code{"-P"} and \code{"-PA"} work in opposite direction of \code{"P"} and \code{"PA"}, respectively.
#'
#' @param msel : If \code{"R"} (default), random selection on males;
#'  if \code{"P"}, selection on phenotypes or true breeding values if \code{Ve} = 0;
#' if \code{"PA"}, selection on true parent averages.
#' \code{"-P"} and \code{"-PA"} work in opposite direction of \code{"P"} and \code{"PA"}, respectively.
#'
#' @return ped : The output pedigree \code{data.frame}. Further information provided in \strong{Details}.
#'
#' @details
#' The output pedigree \code{data.frame} (\code{ped}) has 9 columns: ID, SIRE, DAM, SEX,
#' GEN (generation number starting with 0 for the base generation), PA (parent average),
#' MS (Mendelian Sampling), E (environment and residuals), and P (phenotype).
#'
#' Assortative and disassortative matings can be simulated with different combinations of
#' \code{fsel} and \code{msel}, where both parameters are not set to \code{"R"}.
#'
#' @examples
#' ped = simulatePed(
#'     F0size = 100,
#'     Va0 = 9,
#'     Ve = 36,
#'     littersize = 2,
#'     ngen = 4,
#'     mort.rate = 0.05,
#'     overlap.s = 1,
#'     overlap.d = 0,
#'     f.rate = 0.8,
#'     m.rate = 0.5,
#'     fsel = "P",
#'     msel = "PA"
#' )
#'
#' @export
simulatePed <- function(F0size, Va0, Ve, littersize=1, ngen, mort.rate=0, overlap.s=0, overlap.d=0, f.rate=1, m.rate=1, fsel="R", msel="R") {
    # Check inputs
    ## Check F0size
    F0size = round(F0size)
    if(F0size < 2) stop("ERROR: F0size < 2")
    if((F0size %% 2)!=0) stop("ERROR: F0size should be an even number.")
    ## Check Va0
    if(Va0 <= 0) stop("ERROR: Va0 <= 0")
    ## Check Ve
    if(Ve < 0) stop("ERROR: Ve < 0")
    ## Check littersize
    littersize = round(littersize)
    if(littersize < 1) stop("ERROR: littersize < 1")
    ## Check ngen
    ngen = round(ngen)
    if(ngen < 1) stop("ERROR: ngen < 1")
    ## Check mort.rate
    if(mort.rate > 0.5) stop("ERROR: mort.rate > 0.5")
    if(mort.rate < 0) stop("ERROR: mort.rate < 0")
    ## Check overlap.s
    overlap.s = round(overlap.s)
    if(overlap.s < 0) stop("ERROR: overlap.s < 0")
    if(overlap.s >= ngen) stop("ERROR: overlap.s >= ngen")
    ## Check overlap.d
    overlap.d = round(overlap.d)
    if(overlap.d < 0) stop("ERROR: overlap.d < 0")
    if(overlap.d >= ngen) stop("ERROR: overlap.d >= ngen")
    ## Check f.rate
    if(f.rate > 1) stop("ERROR: f.rate > 1")
    if(f.rate <= 0) stop("ERROR: f.rate <= 0")
    ## Check m.rate
    if(m.rate > 1) stop("ERROR: m.rate > 1")
    if(m.rate <= 0) stop("ERROR: m.rate <= 0")
    if(m.rate > f.rate) stop("ERROR: m.rate > f.rate")
    ## Check fsel
    if(!fsel %in% c("R","P","PA","-P","-PA")) stop('ERROR: fsel should be "R", "P", "PA", "-P" or "-PA".')
    ## Check msel
    if(!msel %in% c("R","P","PA","-P","-PA")) stop('ERROR: msel should be "R", "P", "PA", "-P" or "-PA".')
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
            # Mortality after maturity is applied via overlap.s & overlap.d
            if(mort.rate > 0) {
                ndead = round(mort.rate*nrow(ped[ped$GEN==(i-1),]))
                ped[sample(which(ped$GEN==(i-1)), ndead), "DEAD"] = TRUE
            }
            # Find selection candidates
            sires = ped[ped$SEX=="m" & ped$GEN %in% (-overlap.s:0)+i-1 & !ped$DEAD, c("ID","SBV","DBV","P")]
            dams  = ped[ped$SEX=="f" & ped$GEN %in% (-overlap.d:0)+i-1 & !ped$DEAD, c("ID","SBV","DBV","P")]
            # Selection on males
            nm = round(m.rate*nrow(sires))
            if(nm==0) stop("ERROR: No male left. Generation ", i)
            if(msel=="R") {
                sires = sires[sample(1:nrow(sires), nm),]$ID
            } else if(msel=="P") {
                sires = sires[order(-sires$P),]$ID[1:nm]
            } else if(msel=="PA") {
                sires = sires[order(-sires$SBV -sires$DBV),]$ID[1:nm]
            } else if(msel=="-P") {
                sires = sires[order(sires$P),]$ID[1:nm]
            } else if(msel=="-PA") {
                sires = sires[order(sires$SBV + sires$DBV),]$ID[1:nm]
            }
            # Selection on females
            nf = round(f.rate*nrow(dams))
            if(fsel=="R") {
                dams = dams[sample(1:nrow(dams), nf),]$ID
            } else if(fsel=="P") {
                dams = dams[order(-dams$P),]$ID[1:nf]
            } else if(fsel=="PA") {
                dams = dams[order(-dams$SBV -dams$DBV),]$ID[1:nf]
            } else if(fsel=="-P") {
                dams = dams[order(dams$P),]$ID[1:nf]
            } else if(fsel=="-PA") {
                dams = dams[order(dams$SBV + dams$DBV),]$ID[1:nf]
            }
            # Set mates
            SIRES = c(rep(sires, each=floor(nf/nm)), sires[0:(nf-length(sires)*floor(nf/nm))])
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
