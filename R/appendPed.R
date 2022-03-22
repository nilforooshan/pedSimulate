#' @title Simulate new generations from an existing pedigree
#'
#' @description Simulate pedigree, genetic merits and phenotypes with random/assortative/disassortative matings
#' followed by random/non-random selection of males and females with similar/different selection patterns in males and females,
#' starting from an existing pedigree.
#'
#' @param ped : The input pedigree \code{data.frame} with 9 columns: ID, SIRE, DAM, SEX,
#' GEN (generation), PA (parent average), MS (Mendelian Sampling), E (environment and residuals), and P (phenotype).
#' Note that PA + MS + E = P - \eqn{\mu}, where \eqn{\mu} is the population mean,
#' and PA + MS = BV (genetic merit or breeding value).
#' If MS and E are unknown, those can be set to 0.
#' PA should be equal to the average of sire BV (SBV) and dam BV (DBV).
#' If this condition is not met, PA - (SBV + DBV)/2 is added to MS and (SBV + DBV)/2 replaces PA.
#'
#' @param Va0 : Additive genetic variance in the base generation (i.e., F0).
#'
#' @param Ve : Residual variance, constant across generations.
#'
#' @param littersize : Litter size, default = 1.
#'
#' @param ngen : Number of generations to simulate after the founder generation.
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
#' @param f.order : Ordering selected females for mating;
#' if \code{"fsel"} (default), same as the selection order;
#' if \code{"R"} random ordering;
#' if \code{"P"}, ordering based on phenotypes or true breeding values if \code{Ve} = 0;
#' if \code{"PA"}, ordering based on true parent averages.
#' \code{"-P"} and \code{"-PA"} work in opposite direction of \code{"P"} and \code{"PA"}, respectively.
#'
#' @param m.order : Ordering selected males for mating;
#' if \code{"msel"} (default), same as the selection order;
#' if \code{"R"} random ordering;
#' if \code{"P"}, ordering based on phenotypes or true breeding values if \code{Ve} = 0;
#' if \code{"PA"}, ordering based on true parent averages.
#' \code{"-P"} and \code{"-PA"} work in opposite direction of \code{"P"} and \code{"PA"}, respectively.
#'
#' @param seed : A numeric variable input to the random number generator for reproducible simulations,
#' default = `NA` for non-reproducible simulations.
#'
#' @return ped2 : New generations appended to the input pedigree \code{data.frame}.
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
#' ped2 = appendPed(
#'     ped = ped,
#'     Va0 = 9,
#'     Ve = 36,
#'     littersize = 2,
#'     ngen = 2,
#'     mort.rate = 0.05,
#'     overlap.s = 1,
#'     overlap.d = 0,
#'     f.rate = 0.8,
#'     m.rate = 0.5,
#'     fsel = "R",
#'     msel = "R",
#'     f.order = "P",
#'     m.order = "PA",
#'     seed = 76
#' )
#'
#' @export
appendPed <- function(ped, Va0, Ve, littersize=1, ngen, mort.rate=0, overlap.s=0, overlap.d=0, f.rate=1, m.rate=1, fsel="R", msel="R", f.order="fsel", m.order="msel", seed=NA) {
    # Check inputs
    ## Check ped
    if(!identical(colnames(ped), c("ID","SIRE","DAM","SEX","GEN","PA","MS","E","P"))) stop('colnames(ped) is not c("ID","SIRE","DAM","SEX","GEN","PA","MS","E","P")')
    if(!identical(as.integer(ped[,1]), 1:nrow(ped))) stop("!identical(as.integer(ped[,1]), 1:nrow(ped)); Please consider renumberring the pedigree: ped[,1:3] <- ggroups::renum(ped[,1:3])$newped")
    ## Check Va0
    stopifnot(Va0 > 0)
    ## Check Ve
    stopifnot(Ve >= 0)
    ## Check littersize
    littersize = round(littersize)
    stopifnot(littersize >= 1)
    ## Check ngen
    ngen = round(ngen)
    stopifnot(ngen >= 1)
    ## Check mort.rate
    stopifnot(mort.rate <= 0.5)
    stopifnot(mort.rate >= 0)
    ## Check overlap.s
    overlap.s = round(overlap.s)
    stopifnot(overlap.s >= 0)
    ## Check overlap.d
    overlap.d = round(overlap.d)
    stopifnot(overlap.d >= 0)
    ## Check f.rate
    stopifnot(f.rate <= 1)
    stopifnot(f.rate > 0)
    ## Check m.rate
    stopifnot(m.rate <= 1)
    stopifnot(m.rate > 0)
    stopifnot(m.rate <= f.rate)
    ## Check fsel
    stopifnot(fsel %in% c("R","P","PA","-P","-PA"))
    ## Check msel
    stopifnot(msel %in% c("R","P","PA","-P","-PA"))
    ## Check f.order
    stopifnot(f.order %in% c("fsel", "R","P","PA","-P","-PA"))
    if(f.order=="fsel") f.order = fsel
    ## Check m.order
    stopifnot(m.order %in% c("msel", "R","P","PA","-P","-PA"))
    if(m.order=="msel") m.order = msel
    stopifnot(length(seed)==1)
    stopifnot(is.na(seed) | is.numeric(seed))
    if(!is.na(seed)) set.seed(seed)
    ## Report what you got
    message("ped has ", nrow(ped), " rows.\n",
            "Va0 = ", Va0, "\n",
            "Ve = ", Ve, "\n",
            "littersize = ", littersize, "\n",
            "ngen = ", ngen, "\n",
            "mort.rate = ", mort.rate, "\n",
            "overlap.s = ", overlap.s, "\n",
            "overlap.d = ", overlap.d, "\n",
            "f.rate = ", f.rate, "\n",
            "m.rate = ", m.rate, "\n",
            "fsel = ", fsel, "\n",
            "msel = ", msel, "\n",
            "f.order = ", f.order, "\n",
            "m.order = ", m.order, "\n",
            "seed = ", seed)
    ordcol = c("ID","SIRE","DAM","SEX","GEN","SBV","DBV","MS","E","P","DEAD")
    # Process the input pedigree
    ped$BV = ped$PA + ped$MS
    ped = merge(ped, ped[,c("ID","BV")], by.x="SIRE", by.y="ID", all.x=TRUE)
    colnames(ped)[which(colnames(ped)=="BV.x")] = "BV"
    colnames(ped)[which(colnames(ped)=="BV.y")] = "SBV"
    ped = merge(ped, ped[,c("ID","BV")], by.x="DAM", by.y="ID", all.x=TRUE)
    colnames(ped)[which(colnames(ped)=="BV.x")] = "BV"
    colnames(ped)[which(colnames(ped)=="BV.y")] = "DBV"
    ped[is.na(ped$SBV),]$SBV = 0
    ped[is.na(ped$DBV),]$DBV = 0
    ped$MS = ped$MS + ped$PA - (ped$SBV + ped$DBV)/2
    ped$DEAD = FALSE
    ped = ped[order(ped$ID),]
    maxGEN = max(ped$GEN)
    minGEN = min(ped$GEN)
    message("Started with a pedigree of ", nrow(ped), " individuals and maximum generation number ", maxGEN)
    # Find the population mean
    popmean = mean(ped$P - ped$PA - ped$MS - ped$E)
    ped = ped[,ordcol]
    # Mortality before the last generation
    if(mort.rate > 0 & minGEN!=maxGEN) {
        for(i in (minGEN+1):maxGEN)
        {
            ndead = round(mort.rate*nrow(ped[ped$GEN==(i-1),]))
            if(ndead > 0) ped[sample(which(ped$GEN==(i-1)), ndead), "DEAD"] = TRUE
        }
    }
    # Simulating the next generations
    for(i in maxGEN+(1:ngen))
    {
        # Mortality before maturity
        # Mortality after maturity is applied via overlap.s & overlap.d
        if(mort.rate > 0) {
            ndead = round(mort.rate*nrow(ped[ped$GEN==(i-1),]))
            if(ndead > 0) ped[sample(which(ped$GEN==(i-1)), ndead), "DEAD"] = TRUE
        }
        # Find selection candidates
        sires = ped[ped$SEX=="m" & ped$GEN %in% ((-overlap.s:0)+i-1) & !ped$DEAD, c("ID","SBV","DBV","P")]
        dams  = ped[ped$SEX=="f" & ped$GEN %in% ((-overlap.d:0)+i-1) & !ped$DEAD, c("ID","SBV","DBV","P")]
        # Selection on males
        nm = round(m.rate*nrow(sires))
        if(nm==0) stop("No male left. Generation ", i)
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
        # Order males
        if(m.order!=msel) {
            sires = ped[ped$ID %in% sires,]
            if(m.order=="R") {
                if(nrow(sires) > 1) {
                    sires = sample(sires$ID)
                } else {
                    sires = sires$ID
                }
            } else if(m.order=="P") {
                sires = sires[order(-sires$P),]$ID
            } else if(m.order=="PA") {
                sires = sires[order(-sires$SBV -sires$DBV),]$ID
            } else if(m.order=="-P") {
                sires = sires[order(sires$P),]$ID
            } else if(m.order=="-PA") {
                sires = sires[order(sires$SBV + sires$DBV),]$ID
            }
        }
        # Order females
        if(f.order!=fsel) {
            dams = ped[ped$ID %in% dams,]
            if(f.order=="R") {
                if(nrow(dams) > 1) {
                    dams = sample(dams$ID)
                } else {
                    dams = dams$ID
                }
            } else if(f.order=="P") {
                dams = dams[order(-dams$P),]$ID
            } else if(f.order=="PA") {
                dams = dams[order(-dams$SBV -dams$DBV),]$ID
            } else if(f.order=="-P") {
                dams = dams[order(dams$P),]$ID
            } else if(f.order=="-PA") {
                dams = dams[order(dams$SBV + dams$DBV),]$ID
            }
        }
        # Set mates
        SIRES = rep(sires, each=floor(nf/nm))
        if(length(SIRES) < length(dams)) SIRES = c(SIRES, sires[1:(nf-length(sires)*floor(nf/nm))])
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
        tmp$P = (tmp$SBV + tmp$DBV)/2 + tmp$MS + tmp$E + popmean
        tmp$DEAD = FALSE
        message("Simulated generation ", i, " with ", nrow(tmp), " individuals.")
        ped = rbind(ped, tmp[,ordcol])
    }
    # Replace SBV & DBV with PA
    ped$PA = (ped$SBV + ped$DBV)/2
    ped = ped[,c("ID","SIRE","DAM","SEX","GEN","PA","MS","E","P")]
    ped = ped[order(ped$ID),]
    rownames(ped) = 1:nrow(ped)
    return(ped)
}
