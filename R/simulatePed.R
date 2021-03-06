#' @title Simulate pedigree, genetic merits and phenotypes
#'
#' @description Simulate Pedigree, genetic merits and phenotypes with random mating
#' followed by (non)random selection with different patterns in males and females.
#'
#' @param F0size : Even number of founder animals. No mortality and selection in this generation.
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
#' @param fullsib : If \code{FALSE}, avoid fullsib matings, default = \code{TRUE}.
#' Further information provided in \strong{Details}.
#'
#' @param parentprogeny : If \code{FALSE}, avoid parent-progeny matings, default = \code{TRUE}.
#' Further information provided in \strong{Details}.
#'
#' @param Va0 : Additive genetic variance in the base generation (i.e., F0).
#'
#' @param Ve : Environment (plus residual) variance, set constant across generations.
#'
#' @param fsel : If \code{"R"} (default), random selection on females;
#' if \code{"P"}, selection on phenotypes or true breeding values if \code{Ve} = 0;
#' if \code{"PA"}, selection on true parent averages; redundant if \code{f.rate = 1}.
#' \code{"-P"} and \code{"-PA"} work in opposite direction of \code{"P"} and \code{"PA"}, respectively.
#'
#' @param msel : If \code{"R"} (default), random selection on males;
#'  if \code{"P"}, selection on phenotypes or true breeding values if \code{Ve} = 0;
#' if \code{"PA"}, selection on true parent averages; redundant if \code{m.rate = 1}.
#' \code{"-P"} and \code{"-PA"} work in opposite direction of \code{"P"} and \code{"PA"}, respectively.
#'
#' @return ped : The output pedigree \code{data.frame}. Further information provided in \strong{Details}.
#'
#' @details
#' \code{fullsib = FALSE} : Avoid fullsib matings in each generation by replacing the male mate (SIRE)
#' with a random SIRE among the selected sires, until no fullsib mating is left.
#' If due to a small population bottleneck there is any fullsib mating remained, it would be reported.
#'
#' \code{parentprogeny = FALSE} : Avoid parent-progeny matings in each generation by replacing the male mate (SIRE)
#' with a random SIRE among the selected sires, until no parent-progeny mating is left.
#' If due to a small population bottleneck there is any parent-progeny mating remained, it would be reported.
#'
#' The output pedigree \code{data.frame} (\code{ped}) has 9 columns: ID, SIRE, DAM, SEX,
#' GEN (generation number starting with 0 for the base generation), PA (parent average),
#' MS (Mendelian Sampling), E (environment and residuals), and P (phenotype).
#'
#' @examples
#' ped = simulatePed(
#'     F0size = 100,
#'     f.rate = 0.8,
#'     m.rate = 0.5,
#'     mort.rate = 0.05,
#'     littersize = 2,
#'     ngen = 4,
#'     overlap.s = 1,
#'     overlap.d = 0,
#'     fullsib = FALSE,
#'     parentprogeny = FALSE,
#'     Va0 = 9,
#'     Ve = 36,
#'     fsel = "P",
#'     msel = "PA"
#' )
#'
#' @export
simulatePed <- function(F0size, f.rate=1, m.rate=1, mort.rate=0, littersize=1, ngen, overlap.s=0, overlap.d=0, fullsib=TRUE, parentprogeny=TRUE, Va0, Ve, fsel="R", msel="R") {
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
    if(!fullsib %in% c(TRUE, FALSE)) stop("ERROR: fullsib should be TRUE or FALSE.")
    if(littersize < 1) stop("ERROR: littersize < 1")
    if(ngen < 1) stop("ERROR: ngen < 1")
    if(overlap.s < 0) stop("ERROR: overlap.s < 0")
    if(overlap.d < 0) stop("ERROR: overlap.d < 0")
    if(overlap.s >= ngen) stop("ERROR: overlap.s >= ngen")
    if(overlap.d >= ngen) stop("ERROR: overlap.d >= ngen")
    if(!parentprogeny %in% c(TRUE, FALSE)) stop("ERROR: parentprogeny should be TRUE or FALSE.")
    if(overlap.s==0 & overlap.d==0) parentprogeny = TRUE
    if(Va0 < 0) stop("ERROR: Va0 < 0")
    if(Ve < 0) stop("ERROR: Ve < 0")
    if(!fsel %in% c("R","P","PA","-P","-PA")) stop('ERROR: fsel should be "R", "P", "PA", "-P" or "-PA".')
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
            if(mort.rate > 0) {
                ndead = round(mort.rate*nrow(ped[ped$GEN==(i-1),]))
                ped[sample(which(ped$GEN==(i-1)), ndead), "DEAD"] = TRUE
            }
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
            tmp = data.frame(SIRE=c(rep(sires, floor(nf/nm)), sample(sires, nf-length(sires)*floor(nf/nm))),
                             DAM=sample(dams, length(dams)))
            # Avoid fullsib matings
            if(!fullsib) {
                fs.mate = fs_mate_finder(tmp, ped[,1:3])
                if(nrow(fs.mate) > 0) {
                    prev.fs = nrow(fs.mate) + 1
                    while(nrow(fs.mate) < prev.fs)
                    {
                        prev.fs = nrow(fs.mate)
                        fs.mates = paste(tmp$SIRE, tmp$DAM) %in% paste(fs.mate[,1], fs.mate[,2])
                        tmp[fs.mates,]$SIRE = sample(sires, nrow(tmp[fs.mates,]))
                        fs.mate = fs_mate_finder(tmp[fs.mates,], ped[,1:3])
                    }
                }
            }
            # Avoid parent-progeny matings
            if(!parentprogeny) {
                pp.mate = pp_mate_finder(tmp, ped[,1:3])
                if(nrow(pp.mate) > 0) {
                    prev.pp = nrow(pp.mate) + 1
                    while(nrow(pp.mate) < prev.pp)
                    {
                        prev.pp = nrow(pp.mate)
                        pp.mates = paste(tmp$SIRE, tmp$DAM) %in% paste(pp.mate[,1], pp.mate[,2])
                        tmp[pp.mates,]$SIRE = sample(sires, nrow(tmp[pp.mates,]))
                        pp.mate = fs_mate_finder(tmp[pp.mates,], ped[,1:3])
                    }
                }
            }
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
    # Report fullsib matings, if any
    if(!fullsib) {
        fs.mate = fs_mate_finder(unique(ped[,2:3]), ped[,1:3])
        if(nrow(fs.mate) > 0) {
            warning("Found fullsib matings.")
            for(j in 1:nrow(fs.mate)) message(fs.mate[j,1], " ", fs.mate[j,2])
        }
    }
    # Report parent-progeny matings, if any
    if(!parentprogeny) {
        pp.mate = fs_mate_finder(unique(ped[,2:3]), ped[,1:3])
        if(nrow(pp.mate) > 0) {
            warning("Found parent-progeny matings.")
            for(j in 1:nrow(pp.mate)) message(pp.mate[j,1], " ", pp.mate[j,2])
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
