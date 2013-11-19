################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Methods for objects of class "twinstim", specifically:
### vcov, logLik, print, summary, plot (intensity, iaf), R0, residuals, update
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 666 $
### $Date: 2013-11-08 15:45:36 +0100 (Fre, 08 Nov 2013) $
################################################################################


### don't need a specific coef-method (identical to stats:::coef.default)
## coef.twinstim <- function (object, ...)
## {
##     object$coefficients
## }

### but a method to extract the coefficients in a list might be useful
## "npars" is a named vector or list, where each element states the number of
## coefficients in a group of coefficients.
coeflist <- function (coefs, npars)
{
    coeflist <- vector("list", length(npars))
    names(coeflist) <- names(npars)
    csum <- c(0L,cumsum(npars))
    for (i in seq_along(coeflist)) {
        coeflist[[i]] <- coefs[seq.int(from=csum[i]+1L, length.out=npars[[i]])]
    }
    coeflist
}

## extended and twinstim-specific version,
## which renames elements and unions nbeta0 and p as "endemic"
mycoeflist <- function (coefs, npars)
{
    coeflist <- coeflist(coefs, npars)
    coeflist <- c(list(c(coeflist[[1]], coeflist[[2]])), coeflist[-(1:2)])
    names(coeflist) <- c("endemic", "epidemic", "siaf", "tiaf")
    coeflist
}


## asymptotic variance-covariance matrix (inverse of fisher information matrix)
vcov.twinstim <- function (object, ...)
{
    solve(object$fisherinfo)  # inverse of estimated expected fisher information
}

## Extract log-likelihood of the model (which also enables the use of AIC())
logLik.twinstim <- function (object, ...)
{
    r <- object$loglik
    attr(r, "df") <- length(coef(object))
    class(r) <- "logLik"
    r
}

## Also define an extractAIC-method to make step() work
extractAIC.twinstim <- function (fit, scale, k = 2, ...)
{
    loglik <- logLik(fit)
    edf <- attr(loglik, "df")
    penalty <- k * edf
    c(edf = edf, AIC = -2 * c(loglik) + penalty)            
}

## Number of events (excluding the pre-history)
nobs.twinstim <- function (object, ...) length(object$fitted)

## print-method
print.twinstim <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")
    print.default(x$call)
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits=digits), print.gap = 2, quote = FALSE)
    cat("\nLog-likelihood: ", format(logLik(x), digits=digits), "\n", sep = "")
    if (!isTRUE(x$converged)) {
        cat("\nWARNING: OPTIMIZATION ROUTINE DID NOT CONVERGE!",
            paste0("(",x$converged,")"), "\n")
    }
    cat("\n")
    invisible(x)
}

summary.twinstim <- function (object, test.iaf = FALSE,
    correlation = FALSE, symbolic.cor = FALSE, ...)
{
    ans <- unclass(object)[c("call", "converged", "counts")]
    ans$cov <- vcov(object)
    npars <- object$npars
    coefs <- coef(object)
    nbeta0 <- npars[1]; p <- npars[2]; nbeta <- nbeta0 + p
    q <- npars[3]
    nNotIaf <- nbeta + q
    niafpars <- npars[4] + npars[5]
    est <- coefs
    se <- sqrt(diag(ans$cov))
    zval <- est/se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    coefficients <- cbind(est, se, zval, pval)
    dimnames(coefficients) <- list(names(est),
        c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    ans$coefficients.beta <- coefficients[seq_len(nbeta),,drop=FALSE]
    ans$coefficients.gamma <- coefficients[nbeta+seq_len(q),,drop=FALSE]
    ans$coefficients.iaf <- coefficients[nNotIaf+seq_len(niafpars),,drop=FALSE]
    if (!test.iaf) {
        ## usually, siaf and tiaf parameters are strictly positive,
        ## or parametrized on the logscale. In this case the usual wald test
        ## with H0: para=0 is invalid or meaningless.
        is.na(ans$coefficients.iaf[,3:4]) <- TRUE
    }
    # estimated parameter correlation
    if (correlation) {
        ans$correlation <- cov2cor(ans$cov)
        ans$symbolic.cor <- symbolic.cor
    }
    ans$loglik <- logLik(object)
    ans$aic <- AIC(object)
    ans$runtime <- object$runtime
    class(ans) <- "summary.twinstim"
    ans
}

## additional methods to make confint.default work for summary.twinstim
vcov.summary.twinstim <- function (object, ...) object$cov
coef.summary.twinstim <- function (object, ...) with(object, {
    coeftab <- rbind(coefficients.beta, coefficients.gamma, coefficients.iaf)
    structure(coeftab[,1], names=rownames(coeftab))
})

## print-method for summary.twinstim
print.summary.twinstim <- function (x,
    digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...)
{
    nbeta <- nrow(x$coefficients.beta) # = nbeta0 + p
    q <- nrow(x$coefficients.gamma)
    niafpars <- nrow(x$coefficients.iaf)
    cat("\nCall:\n")
    print.default(x$call)
    if (nbeta > 0L) {
        cat("\nCoefficients of the endemic component:\n")
        printCoefmat(x$coefficients.beta, digits = digits,
            signif.stars = signif.stars, signif.legend = (q==0L) && signif.stars, ...)
    } else cat("\nNo coefficients in the endemic component.\n")
    if (q + niafpars > 0L) {
        cat("\nCoefficients of the epidemic component:\n")
        printCoefmat(rbind(x$coefficients.gamma, x$coefficients.iaf), digits = digits,
            signif.stars = signif.stars, ...)
#         if (niafpars > 0L) {
#             #cat("Coefficients of interaction functions:\n")
#             printCoefmat(x$coefficients.iaf, digits = digits, signif.stars = signif.stars, ...)
#         }
    } else cat("\nNo epidemic component.\n")
    cat("\nAIC: ", format(x$aic, digits=max(4, digits+1)))
    cat("\nLog-likelihood:", format(x$loglik, digits = digits))
    cat("\nNumber of log-likelihood evaluations:", x$counts[1])
    cat("\nNumber of score function evaluations:", x$counts[2])
    cores <- attr(x$runtime, "cores")
    cat("\nRuntime",
        if (!is.null(cores) && cores > 1) paste0(" (", cores, " cores)")
        , ": ", format(
        if (length(x$runtime)==0) NA_real_ else if (length(x$runtime)==1)
        x$runtime else x$runtime[["elapsed"]],
        digits=max(4, digits+1)), " seconds", sep="")
    cat("\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1L) {
        cat("\nCorrelation of Coefficients:\n")
        if (is.logical(symbolic.cor) && symbolic.cor) {
            correl <- symnum(correl, abbr.colnames = NULL)
            correlcodes <- attr(correl, "legend")
            attr(correl, "legend") <- NULL
            print(correl)
            cat("---\nCorr. codes:  ", correlcodes, "\n", sep="")
        } else {
            correl <- format(round(correl, 2), nsmall = 2)
            correl[!lower.tri(correl)] <- ""
            colnames(correl) <- substr(colnames(correl), 1, 5)
            print(correl[-1, -p, drop = FALSE], quote = FALSE)
        }
        }
    }
    if (!isTRUE(x$converged)) {
        cat("\nWARNING: OPTIMIZATION ROUTINE DID NOT CONVERGE!",
            paste0("(",x$converged,")"), "\n")
    }
    cat("\n")
    invisible(x)
}



### 'cat's the summary in LaTeX code

toLatex.summary.twinstim <- function (object,
                                      digits = max(3, getOption("digits") - 3),
                                      eps.Pvalue = 1e-4,
                                      align = "lrrrr", withAIC = FALSE, ...)
{
ret <- capture.output({
    cat("\\begin{tabular}{", align, "}\n\\hline\n", sep="")
    cat(" & Estimate & Std. Error & $z$ value & $P(|Z|>|z|)$ \\\\\n\\hline\n\\hline\n")
    tabh <- object$coefficients.beta
    tabe <- rbind(object$coefficients.gamma, object$coefficients.iaf)
    for (tabname in c("tabh", "tabe")) {
        tab <- get(tabname)
        if (nrow(tab) > 0L) {
            rownames(tab) <- gsub(" ", "", rownames(tab))
            tab_char <- capture.output(
                        printCoefmat(tab, digits=digits, signif.stars=FALSE,
                                     eps.Pvalue = eps.Pvalue, na.print="NA")
                        )[-1]
            #tab_char <- sub("< (0\\..+)$", "<\\1", tab_char)
            tab_char <- sub("([<]?)[ ]?([0-9]+)e([+-][0-9]+)$",
                            "\\1\\2\\\\cdot{}10^{\\3}",
                            tab_char)
            con <- textConnection(tab_char)
            tab2 <- read.table(con, colClasses="character")
            close(con)
            rownames(tab2) <- paste0("\\texttt{",gsub("_","\\\\_",tab2[,1]),"}")
            tab2 <- tab2[,-1]
            tab2[] <- lapply(tab2, function(x) {
                ifelse(is.na(x), "", paste0("$",x,"$")) # (test.iaf=FALSE)
            })
            print(xtable(tab2), only.contents=TRUE, hline.after=NULL,
                  include.colnames=FALSE, sanitize.text.function=identity)
            cat("\\hline\n")
        }
    }
    if (withAIC) {
        cat("\\hline\n")
        cat("AIC:& $", format(object$aic, digits=max(4, digits+1)), "$ &&&\\\\\n")
        cat("Log-likelihood:& $", format(object$loglik, digits=digits), "$ &&&\\\\\n")
    }
    cat("\\hline\n")
    cat("\\end{tabular}\n")
})
class(ret) <- "Latex"
ret
}


## Alternative implementation including exp-transformation of parameters
## CAVE: no iaf parameters here

xtable.summary.twinstim <- function (x, caption = NULL, label = NULL,
                             align = c("l", "r", "r", "r"), digits = 3,
                             display = c("s", "f", "s", "s"),
                             ci.level = 0.95, ci.fmt = "%4.2f", ci.to = "--",
                             eps.pvalue = 1e-4, ...)
{
    cis <- confint(x, level=ci.level)
    tabh <- x$coefficients.beta
    tabe <- x$coefficients.gamma
    tab <- rbind(tabh, tabe)
    tab <- tab[grep("^([he]\\.\\(Intercept\\)|h.type)", rownames(tab),
                    invert=TRUE),,drop=FALSE]
    expcis <- exp(cis[rownames(tab),,drop=FALSE])
    cifmt <- paste0(ci.fmt, ci.to, ci.fmt)
    rrtab <- data.frame(RR = exp(tab[,1]),
                        CI = sprintf(cifmt, expcis[,1], expcis[,2]),
                        "p-value" = formatPval(tab[,4], eps=eps.pvalue),
                        check.names = FALSE, stringsAsFactors=FALSE)
    names(rrtab)[2] <- paste0(100*ci.level, "% CI")

    ## append caption etc.
    class(rrtab) <- c("xtable", "data.frame")
    caption(rrtab) <- caption
    label(rrtab) <- label
    align(rrtab) <- align
    digits(rrtab) <- digits
    display(rrtab) <- display

    ## Done
    rrtab
}

xtable.twinstim <- function () {
    cl <- match.call()
    cl[[1]] <- as.name("xtable.summary.twinstim")
    cl$x <- substitute(summary(x))
    eval.parent(cl)                # => xtable.summary.twinstim must be exported
}
formals(xtable.twinstim) <- formals(xtable.summary.twinstim)


### Plot temporal or spatial evolution of the intensity

intensity.twinstim <- function (x, aggregate = c("time", "space"),
    types = 1:nrow(x$qmatrix), tiles, tiles.idcol = NULL)
{
    modelenv <- environment(x)
    
    ## check arguments
    if (is.null(modelenv))
        stop("'x' is missing the model environment\n",
             "  -- re-fit or update() with 'model=TRUE'")
    aggregate <- match.arg(aggregate)
    stopifnot(is.vector(types, mode="numeric"),
              types %in% seq_len(modelenv$nTypes),
              !anyDuplicated(types))

    ## remove (big) x object from current evaluation environment
    qmatrix <- x$qmatrix                # not part of modelenv
    force(types)                        # evaluate types before rm(x)
    rm(x)                               # don't need this anymore

    ##thisenv <- environment()
    ##parent.env(thisenv) <- modelenv     # objects of modelenv become visible
    ## Instead of the above, we do cheap and nasty model unpacking!
    ## safer than the parent.env<- hack (R manual: "extremely dangerous"), and
    ## cleaner than running code inside with(modelenv,...) since assignments
    ## then would take place in modelenv, which would produce garbage
    t0 <- modelenv$t0
    T <- modelenv$T
    histIntervals <- modelenv$histIntervals
    eventTimes <- modelenv$eventTimes
    eventCoords <- modelenv$eventCoords
    eventTypes <- modelenv$eventTypes
    removalTimes <- modelenv$removalTimes
    gridTiles <- modelenv$gridTiles
    gridBlocks <- modelenv$gridBlocks
    tiaf <- modelenv$tiaf
    tiafpars <- modelenv$tiafpars
    eps.s <- modelenv$eps.s
    siaf <- modelenv$siaf
    siafpars <- modelenv$siafpars
    
    ## endemic component on the spatial or temporal grid
    hInt <- 
        if (modelenv$hash) {
            eta <- drop(modelenv$mmhGrid %*% modelenv$beta)
            if (!is.null(modelenv$offsetGrid)) eta <- modelenv$offsetGrid + eta
            expeta <- exp(unname(eta))
            .beta0 <- rep_len(if (modelenv$nbeta0==0L) 0 else modelenv$beta0,
                              modelenv$nTypes)
            fact <- sum(exp(.beta0[types]))
            if (aggregate == "time") {      # int over W and types by BLOCK
                fact * c(tapply(expeta * modelenv$ds, gridBlocks, sum,
                                simplify = TRUE))
            } else {                        # int over T and types by tile
                fact * c(tapply(expeta * modelenv$dt, gridTiles, sum,
                                simplify = TRUE))
            }
        } else {
            ngrid <- if (aggregate == "time") {
                gridBlocks[length(gridBlocks)]
            } else nlevels(gridTiles)
            rep.int(0, ngrid)
        }

    ## endemic component as a function of time or location
    hIntFUN <- if (modelenv$hash) {
        if (aggregate == "time") {
            function (tp) {
                stopifnot(isScalar(tp))
                if (tp == t0) hInt[1L] else {
                    starts <- histIntervals$start
                    idx <- match(TRUE, c(starts,T) >= tp) - 1L
                    block <- histIntervals$BLOCK[idx]
                    hInt[as.character(block)]
                }
            }
        } else {
            stopifnot(is(tiles, "SpatialPolygons"))
            tilesIDs <- if (is.null(tiles.idcol)) {
                sapply(tiles@polygons, slot, "ID")
            } else tiles@data[[tiles.idcol]]
            if (!all(levels(gridTiles) %in% tilesIDs)) {
                stop("'tiles' is incomplete for 'x' (check 'tiles.idcol')")
            }
            .tiles <- as(tiles, "SpatialPolygons") # for over-method (drop data)
            function (xy) {             # works with a whole coordinate matrix
                points <- SpatialPoints(xy, proj4string=tiles@proj4string)
                polygonidxOfPoints <- over(points, .tiles)
                points.outside <- is.na(polygonidxOfPoints)
                polygonidxOfPoints[points.outside] <- 1L   # tiny hack
                tilesOfPoints <- if (is.null(tiles.idcol)) {
                    sapply(tiles@polygons[polygonidxOfPoints], slot, "ID")
                } else tiles@data[polygonidxOfPoints,tiles.idcol]
                is.na(tilesOfPoints) <- points.outside     # resolve hack
                hInt[tilesOfPoints]       # index by name
            }
        }
    } else function (...) 0

    ## epidemic component
    eInt <- if (modelenv$hase) {
        qSum_types <- rowSums(qmatrix[,types,drop=FALSE])[eventTypes]
        fact <- qSum_types * modelenv$gammapred
        if (aggregate == "time") {  # as a function of time (int over W & types)
            factS <- fact * modelenv$siafInt
            function (tp) {
                stopifnot(isScalar(tp))
                tdiff <- tp - eventTimes
                infectivity <- qSum_types > 0 & (tdiff > 0) & (removalTimes >= tp)
                if (any(infectivity)) {
                    gsources <- tiaf$g(tdiff[infectivity],
                                       tiafpars,
                                       eventTypes[infectivity])
                    intWj <- factS[infectivity] * gsources
                    sum(intWj)
                } else 0
            }
        } else {  # as a function of location (int over time and types)
            factT <- fact * modelenv$tiafInt
            nEvents <- nrow(eventCoords)
            function (xy) {
                stopifnot(is.vector(xy, mode="numeric"), length(xy) == 2L)
                point <- matrix(xy, nrow=nEvents, ncol=2L, byrow=TRUE)
                sdiff <- point - eventCoords
                proximity <- qSum_types > 0 &
                    .rowSums(sdiff^2, nEvents, 2L) <= eps.s^2
                if (any(proximity)) {
                    fsources <- siaf$f(sdiff[proximity,,drop=FALSE],
                                       siafpars,
                                       eventTypes[proximity])
                    intTj <- factT[proximity] * fsources
                    sum(intTj)
                } else 0
            }
        }
    } else function (...) 0

    ## return component functions
    list(hGrid = hInt, hFUN = hIntFUN, eFUN = eInt,
         aggregate = aggregate, types = types)
}

intensityplot.twinstim <- function (x,
    which = c("epidemic proportion", "endemic proportion", "total intensity"),
    aggregate, types, tiles, tiles.idcol, # arguments of intensity.twinstim;
                                          # defaults are set below
    plot = TRUE, add = FALSE, tgrid = 101, rug.opts = list(),
    sgrid = 128, polygons.args = list(), points.args = list(),
    cex.fun = sqrt, ...)
{
    which <- match.arg(which)
    
    ## set up desired intensities
    cl <- match.call()
    cl <- cl[c(1L, match(names(formals(intensity.twinstim)), names(cl), 0L))]
    cl[[1]] <- as.name("intensity.twinstim")
    components <- eval(cl, envir = parent.frame())
    aggregate <- components$aggregate
    types <- components$types
    
    ## define function to plot
    FUN <- function (tmp) {}
    names(formals(FUN)) <- if (aggregate == "time") "times" else "coords"
    body1 <- if (aggregate == "time") expression(
        hGrid <- sapply(times, components$hFUN, USE.NAMES=FALSE),
        eGrid <- sapply(times, components$eFUN, USE.NAMES=FALSE)
        ) else expression(
            hGrid <- unname(components$hFUN(coords)), # takes whole coord matrix
            eGrid <- apply(coords, 1, components$eFUN)
            )
    body2 <- switch(which,
                    "epidemic proportion" = expression(eGrid / (hGrid + eGrid)),
                    "endemic proportion" = expression(hGrid / (hGrid + eGrid)),
                    "total intensity" = expression(hGrid + eGrid))
    body(FUN) <- as.call(c(as.name("{"), c(body1, body2)))
    
    if (!plot) return(FUN)

    ## plot the FUN
    modelenv <- environment(components$eFUN)$modelenv
    dotargs <- list(...)
    nms <- names(dotargs)
    if (aggregate == "time") {
        ## set up grid of x-values (time points where 'which' will be evaluated)
        tgrid <- if (isScalar(tgrid)) {
            seq(modelenv$t0, modelenv$T, length.out=tgrid)
        } else {
            stopifnot(is.vector(tgrid, mode="numeric"))
            sort(tgrid)
        }
        
        ## calculate 'which' on tgrid
        yvals <- FUN(tgrid)
        
        ## plot it
        if(! "xlab" %in% nms) dotargs$xlab <- "time"
        if(! "ylab" %in% nms) dotargs$ylab <- which
        if(! "type" %in% nms) dotargs$type <- "l"
        if(! "ylim" %in% nms) dotargs$ylim <- {
            if (which == "total intensity") c(0,max(yvals)) else c(0,1)
        }
        do.call(if (add) "lines" else "plot", args=c(alist(x=tgrid, y=yvals), dotargs))
        if (is.list(rug.opts)) {
            if (is.null(rug.opts$ticksize)) rug.opts$ticksize <- 0.02
            if (is.null(rug.opts$quiet)) rug.opts$quiet <- TRUE
            eventTimes.types <- modelenv$eventTimes[modelenv$eventTypes %in% types]
            do.call("rug", args = c(alist(x=eventTimes.types), rug.opts))
        }
        invisible(FUN)
    } else {
        .tiles <- as(tiles, "SpatialPolygons") # we don't need the data here
        
        ## set up grid of coordinates where 'which' will be evaluated
        if (isScalar(sgrid)) {
            sgrid <- maptools::Sobj_SpatialGrid(.tiles, n = sgrid)$SG
            ## ensure that sgrid has exactly the same proj4string as .tiles
            ## since CRS(proj4string(.tiles)) might have modified the string
            sgrid@proj4string <- .tiles@proj4string
        }
        sgrid <- as(sgrid, "SpatialPixels")
        
        ## only select grid points inside W (tiles)
        in.tiles <- !is.na(over(sgrid, .tiles))
        sgrid <- sgrid[in.tiles,]
        
        ## calculate 'which' on sgrid
        yvals <- FUN(coordinates(sgrid))
        sgridy <- SpatialPixelsDataFrame(sgrid, data=data.frame(yvals=yvals),
                                         proj4string=.tiles@proj4string)

        ## define sp.layout
        lobjs <- list()
        if (is.list(polygons.args)) {
            nms.polygons <- names(polygons.args)
            if(! "col" %in% nms.polygons) polygons.args$col <- "darkgrey"
            lobjs <- c(lobjs,
                       list(c(list("sp.polygons", .tiles, first=FALSE),
                              polygons.args)))
        }
        if (is.list(points.args)) {
            eventCoords.types <- modelenv$eventCoords[modelenv$eventTypes %in% types,,drop=FALSE]
            ## eventCoords as Spatial object with duplicates counted and removed
            eventCoords.types <- SpatialPoints(eventCoords.types,
                                               proj4string=.tiles@proj4string,
                                               bbox = .tiles@bbox)
            eventCoords.types <- SpatialPointsDataFrame(eventCoords.types,
                data.frame(mult = multiplicity(eventCoords.types)))
            eventCoords.types <- eventCoords.types[!duplicated(coordinates(eventCoords.types)),]
            points.args <- modifyList(list(pch=1, cex=0.5), points.args)
            pointcex <- cex.fun(eventCoords.types$mult)
            pointcex <- pointcex * points.args$cex
            points.args$cex <- NULL
            lobjs <- c(lobjs,
                       list(c(list("sp.points", eventCoords.types, first=FALSE,
                                   cex=pointcex), points.args)))
        }
        if ("sp.layout" %in% nms) {
            if (!is.list(dotargs$sp.layout[[1]])) { # let sp.layout be a list of lists
                dotargs$sp.layout <- list(dotargs$sp.layout)
            }
            lobjs <- c(lobjs, dotargs$sp.layout)
            dotargs$sp.layout <- NULL
        }

        ## plotit
        if (add) message("'add'ing is not possible with 'aggregate=\"space\"'")
        if (! "xlim" %in% nms) dotargs$xlim <- bbox(.tiles)[1,]
        if (! "ylim" %in% nms) dotargs$ylim <- bbox(.tiles)[2,]
        if (! "scales" %in% nms) dotargs$scales <- list(draw = TRUE)
        do.call("spplot", args=c(alist(sgridy, zcol="yvals", sp.layout=lobjs,
                          checkEmptyRC=FALSE), dotargs))
    }
}

## set default arguments for intensityplot.twinstim from intensity.twinstim
formals(intensityplot.twinstim)[names(formals(intensity.twinstim))] <-
    formals(intensity.twinstim)



### Plot fitted tiaf or siaf(cbind(0, r)), r=distance

iafplot <- function (object, which = c("siaf", "tiaf"),
    types = 1:nrow(object$qmatrix), scaled = FALSE,
    conf.type = if (length(pars) > 1) "bootstrap" else "parbounds",
    conf.level = 0.95, conf.B = 999,
    xgrid = 101, col.estimate = rainbow(length(types)), col.conf = col.estimate,
    alpha.B = 0.15, lwd = c(3,1), lty = c(1,2), xlim = NULL, ylim = NULL,
    add = FALSE, xlab = NULL, ylab = NULL, legend = !add && (length(types) > 1),
    ...)
{
    which <- match.arg(which)
    IAF <- object$formula[[which]][[if (which=="siaf") "f" else "g"]]
    coefs <- coef(object)
    ## epidemic intercept (case of no intercept is not addressed atm)
    if (scaled) {
        idxgamma0 <- match("e.(Intercept)", names(coefs), nomatch=0L)
        if (idxgamma0 == 0L) {
            message("no scaling due to missing epidemic intercept")
            scaled <- FALSE
        }
    } else idxgamma0 <- 0L              # if no scaling, gamma0 is 0-length
    gamma0 <- coefs[idxgamma0]
    ## parameters of the interaction function
    idxiafpars <- grep(paste0("^e\\.",which), names(coefs))
    iafpars <- coefs[idxiafpars]
    ## concatenate parameters
    idxpars <- c(idxgamma0, idxiafpars)
    pars <- c(gamma0, iafpars)
    ## type of confidence band
    force(conf.type)                    # per default depends on 'pars'
    if (length(pars) == 0 || any(is.na(conf.type)) || is.null(conf.type)) {
        conf.type <- "none"
    }
    conf.type <- match.arg(conf.type,
                           choices = c("parbounds", "bootstrap", "none"))

    ## add intercept-factor to IAF if scaled=TRUE
    FUN <- if (scaled) {
        function (x, pars, types) {
            gamma0 <- pars[1L]
            iafpars <- pars[-1L]
            exp(gamma0) * IAF(x, iafpars, types)
        }
    } else IAF
    
    ## grid of x-values (t or ||s||) on which FUN will be evaluated
    if (is.null(xlim)) {
        xmax <- if (add) {
            par("usr")[2] / (if (par("xaxs")=="r") 1.04 else 1)
        } else if (which == "siaf") {
            sqrt(sum((object$bbox[,"max"] - object$bbox[,"min"])^2))
        } else {
            diff(object$timeRange)
        }
        xlim <- c(0, xmax)
    }
    xgrid <- if (isScalar(xgrid)) {
        seq(0, xlim[2], length.out=xgrid)
    } else {
        stopifnot(!is.na(xgrid), is.vector(xgrid, mode="numeric"))
        ## xgrid-specification overrides default xlim
        if (is.null(match.call()$xlim)) xlim <- c(0, max(xgrid))
        sort(xgrid)
    }
    
    ## initialize plotting frame
    if (!add) {
        if (is.null(ylim))
            ylim <- c(0, FUN(if (which=="siaf") cbind(0,0) else 0, pars, 1L))
        if (is.null(xlab)) xlab <- if (which == "siaf") {
            expression("Distance " *
                       x *              # group("||",bold(s)-bold(s)[j],"||") *
                       " from host")
        } else {
            expression("Time " * t * " since infectious")
        }
        if (is.null(ylab)) {
            ylab <- if (which == "siaf") {
                expression(f(x))        # f(group("||",bold(s)-bold(s)[j],"||"))
            } else {
                expression(g(t))
            }
            if (scaled) {
                ylab <- as.expression(as.call(list(quote(paste),
                    quote(e^{gamma[0]}), quote(phantom() %.% phantom()),
                    ylab[[1]])))
            }
        }
        plot(xlim, ylim, type="n", xlab = xlab, ylab = ylab, ...)
    }
    
    for (i in seq_along(types)) {
        ## select parameters on which to evaluate iaf
        parSample <- switch(conf.type, parbounds = {
            cis <- confint(object, idxpars, level=conf.level)
            ## all combinations of parameter bounds
            do.call("expand.grid", as.data.frame(t(cis)))
        }, bootstrap = {
            ## bootstrapping parameter values
            rbind(pars, mvrnorm(conf.B, mu=pars,
                                Sigma=vcov(object)[idxpars,idxpars,drop=FALSE]))
        })
        
        ## add confidence limits
        if (!is.null(parSample)) {
            fvalsSample <- apply(parSample, 1, function (pars)
                                 FUN(if(which=="siaf") cbind(xgrid,0) else xgrid,
                                     pars, types[i]))
            lowerupper <- if (conf.type == "parbounds") {
                t(apply(fvalsSample, 1, range))
            } else { # bootstrapped parameter values
                if (is.na(conf.level)) {
                    stopifnot(alpha.B >= 0, alpha.B <= 1)
                    .col <- col2rgb(col.conf[i], alpha=TRUE)[,1]
                    .col["alpha"] <- round(alpha.B*.col["alpha"])
                    .col <- do.call("rgb", args=c(as.list(.col),
                                           maxColorValue = 255))
                    matlines(x=xgrid, y=fvalsSample, type="l", lty=lty[2],
                             col=.col, lwd=lwd[2]) # returns NULL
                } else {
                    t(apply(fvalsSample, 1, quantile,
                            probs=c(0,conf.level) + (1-conf.level)/2))
                }
            }
            if (!is.null(lowerupper)) {
                matlines(x=xgrid, y=lowerupper, type="l", lty=lty[2],
                         col=col.conf[i], lwd=lwd[2])
            }
        }
        
        ## add point estimate
        lines(x=xgrid,
              y=FUN(if(which=="siaf") cbind(xgrid,0) else xgrid, pars, types[i]),
              lty=lty[1], col=col.estimate[i], lwd=lwd[1])
    }
    
    ## add legend
    if (isTRUE(legend) || is.list(legend)) {
        default.legend <- list(x = "topright",
                               legend = rownames(object$qmatrix)[types],
                               col = col.estimate, lty = lty[1], lwd = lwd[1],
                               bty = "n", cex = 0.9, title="type")
        legend.args <- if (is.list(legend)) {
            modifyList(default.legend, legend)
        } else default.legend
        do.call("legend", legend.args)
    }
    invisible()
}


### Plot method for twinstim (wrapper for iafplot and intensityplot)

plot.twinstim <- function (x, which, ...)
{
    cl <- match.call()
    which <- match.arg(which, choices =
                       c(eval(formals(intensityplot.twinstim)$which),
                         eval(formals(iafplot)$which)))
    FUN <- if (which %in% eval(formals(intensityplot.twinstim)$which))
        "intensityplot" else "iafplot"
    cl[[1]] <- as.name(FUN)
    if (FUN == "iafplot") names(cl)[names(cl) == "x"] <- "object"
    eval(cl, envir = parent.frame())
}



### Calculates the basic reproduction number R0 for individuals
### with marks given in 'newevents'

R0.twinstim <- function (object, newevents, trimmed = TRUE, ...)
{
    ## extract model information
    npars <- object$npars
    if (npars["q"] == 0L) {
        message("no epidemic component in model, returning 0-vector")
        if (missing(newevents)) return(object$R0) else {
            return(structure(rep.int(0, nrow(newevents)),
                             names = rownames(newevents)))
        }
    }
    t0 <- object$timeRange[1L]
    T <- object$timeRange[2L]
    typeNames <- rownames(object$qmatrix)
    nTypes <- length(typeNames)
    types <- seq_len(nTypes)
    form <- formula(object)
    siaf <- form$siaf
    tiaf <- form$tiaf
    coefs <- coef(object)
    tiafpars <- coefs[sum(npars[1:4]) + seq_len(npars["ntiafpars"])]
    siafpars <- coefs[sum(npars[1:3]) + seq_len(npars["nsiafpars"])]
    
    if (missing(newevents)) {
        ## if no newevents are supplied, use original events
        if (trimmed) {                  # already calculated by 'twinstim'
            return(object$R0)
        } else {    # untrimmed version (spatio-temporal integral over R+ x R^2)
            ## extract relevant data from model environment
            if (is.null(modelenv <- environment(object))) {
                stop("need model environment for untrimmed R0 of fitted events\n",
                     " -- re-fit or update() with 'model=TRUE'")
            }
            eventTypes <- modelenv$eventTypes
            eps.t <- modelenv$eps.t
            eps.s <- modelenv$eps.s
            gammapred <- modelenv$gammapred
            names(gammapred) <- names(object$R0) # for names of the result
        }
    } else {
        ## use newevents
        stopifnot(is.data.frame(newevents))
        if (!"time" %in% names(newevents)) {
            stop("missing event \"time\" column in 'newevents'")
        }
        if (any(!c("eps.s", "eps.t") %in% names(newevents))) {
            stop("missing \"eps.s\" or \"eps.t\" columns in 'newevents'")
        }
        stopifnot(is.factor(newevents[["type"]]))
        
        ## subset newevents to timeRange
        .N <- nrow(newevents)
        newevents <- subset(newevents, time + eps.t > t0 & time <= T)
        if (nrow(newevents) < .N) {
            message("subsetted 'newevents' to only include events infectious ",
                    "during 'object$timeRange'")
        }

        ## extract columns
        newevents$type <- factor(newevents[["type"]], levels = typeNames)
        eventTimes <- newevents[["time"]]
        eps.t <- newevents[["eps.t"]]
        eps.s <- newevents[["eps.s"]]
        
        ## calculate gammapred for newevents
        epidemic <- terms(form$epidemic, data = newevents, keep.order = TRUE)
        mfe <- model.frame(epidemic, data = newevents,
                           na.action = na.pass, drop.unused.levels = FALSE)
        mme <- model.matrix(epidemic, mfe)
        gamma <- coefs[sum(npars[1:2]) + seq_len(npars["q"])]
        if (ncol(mme) != length(gamma)) {
            stop("epidemic model matrix has the wrong number of columns ",
                 "(check the variable types in 'newevents' (factors, etc))")
        }
        gammapred <- drop(exp(mme %*% gamma))
        names(gammapred) <- rownames(newevents)

        ## now, convert types of newevents to integer codes
        eventTypes <- as.integer(newevents$type)
    }

    ## qSum
    qSumTypes <- rowSums(object$qmatrix)
    qSum <- unname(qSumTypes[eventTypes])


    ## calculate remaining factors of the R0 formula, i.e. siafInt and tiafInt
    
    if (trimmed) {                      # trimmed R0 for newevents
        
        ## integral of g over the observed infectious periods
        .tiafInt <- .tiafIntFUN()
        gIntUpper <- pmin(T - eventTimes, eps.t)
        gIntLower <- pmax(0, t0 - eventTimes)
        tiafInt <- .tiafInt(tiafpars, from=gIntLower, to=gIntUpper,
                            type=eventTypes, G=tiaf$G)
        ## integral of f over the influenceRegion
        bdist <- newevents[[".bdist"]]
        influenceRegion <- newevents[[".influenceRegion"]]
        if (is.null(influenceRegion)) {
            stop("missing \".influenceRegion\" component in 'newevents'")
        }
        noCircularIR <- if (is.null(bdist)) FALSE else all(eps.s > bdist)
        if (attr(siaf, "constant")) {
            iRareas <- sapply(influenceRegion, area.owin)
            ## will be used by .siafInt()
        } else if (! (is.null(siaf$Fcircle) ||
               (is.null(siaf$effRange) && noCircularIR))) {
            if (is.null(bdist)) {
                stop("missing \".bdist\" component in 'newevents'")
            }
        }
        .siafInt <- .siafIntFUN(siaf, noCircularIR=noCircularIR)
        .siafInt.args <- c(alist(siafpars), object$control.siaf$F)
        siafInt <- do.call(".siafInt", .siafInt.args)
        
    } else {                     # untrimmed R0 for original events or newevents

        if (any(is.infinite(eps.t), is.infinite(eps.s))) {
            message("infinite interaction ranges yield infinite R0 values ",
                    "because 'trimmed = FALSE'")
        }
        
        ## integrals of interaction functions for all combinations of type and
        ## eps.s/eps.t in newevents
        typeTcombis <- expand.grid(type=types, eps.t=unique(eps.t),
                                   KEEP.OUT.ATTRS=FALSE)
        typeTcombis$gInt <-
            with(typeTcombis, tiaf$G(eps.t, tiafpars, type)) -
                tiaf$G(rep.int(0,nTypes), tiafpars, types)[typeTcombis$type]
        
        typeScombis <- expand.grid(type=types, eps.s=unique(eps.s),
                                   KEEP.OUT.ATTRS=FALSE)
        typeScombis$fInt <- apply(typeScombis, MARGIN=1, FUN=function (type_eps.s) {
            type <- type_eps.s[1]
            eps.s <- type_eps.s[2]
            if (is.null(siaf$Fcircle)) {
                do.call(siaf$F, c(list(discpoly(c(0,0), eps.s, class="owin"),
                                       siaf$f, siafpars, type),
                                  object$control.siaf$F))
            } else {
                siaf$Fcircle(eps.s, siafpars, type)
            }
        })
        
        ## match combinations to rows of original events or 'newevents'
        eventscombiidxS <- match(paste(eventTypes,eps.s,sep="."),
                                 with(typeScombis,paste(type,eps.s,sep=".")))
        eventscombiidxT <- match(paste(eventTypes,eps.t,sep="."),
                                 with(typeTcombis,paste(type,eps.t,sep=".")))

        siafInt <- typeScombis$fInt[eventscombiidxS]
        tiafInt <- typeTcombis$gInt[eventscombiidxT]
        
    }

    ## return R0 values
    R0s <- qSum * gammapred * siafInt * tiafInt
    R0s
}



### Extract the "residual process" (cf. Ogata, 1988) of a twinstim, i.e. the
### fitted cumulative intensity of the ground process at the event times.
### "generalized residuals similar to those discussed in Cox and Snell (1968)"

residuals.twinstim <- function (object, ...)
{
  res <- object$tau
  if (is.null(res)) {
      if (is.null(modelenv <- environment(object))) {
          stop("residuals not available; re-fit the model with 'cumCIF = TRUE'")
      } else {
          cat("'", substitute(object), "' was fit with disabled 'cumCIF'",
              " -> calculate it now...\n", sep="")
          res <- with(modelenv, LambdagEvents(cumCIF.pb = TRUE))
          cat("Done.\n")
          try({
              objname <- deparse(substitute(object))
              object$tau <- res
              assign(objname, object, envir = parent.frame())
              cat("Note: added the 'tau' component to '", objname, "' for future use.\n", sep="")
          }, silent = TRUE)
      }
  }
  return(res)
}



######################################################################
# Function to compute estimated and profile likelihood based
# confidence intervals. Heavy computations might be necessary!
#
#Params:
# fitted - output from a fit with twinstim
# profile - list with 4D vector as entries - format:
#               c(index, lower, upper, grid size)
#           where index is the index in the coef vector
#                 lower and upper are the parameter limits (can be NA)
#                 grid size is the grid size of the equally spaced grid
#                 between lower and upper (can be 0)
# alpha - (1-alpha)% profile likelihood CIs are computed.
#         If alpha <= 0 then no CIs are computed
# control - control object to use for optim in the profile loglik computations
#
# Returns:
#  list with profile loglikelihood evaluations on the grid
#  and highest likelihood and wald confidence intervals
######################################################################

profile.twinstim <- function (fitted, profile, alpha = 0.05,
    control = list(fnscale = -1, factr = 1e1, maxit = 100), do.ltildeprofile=FALSE,...)
{
  ## Check that input is ok
  profile <- as.list(profile)
  if (length(profile) == 0L) {
    stop("nothing to do")
  }
  lapply(profile, function(one) {
    if (length(one) != 4L) {
      stop("each profile entry has to be of form ",
           "'c(index, lower, upper, grid size)'")
    }})
  if (is.null(fitted[["functions"]])) {
    stop("'fitted' must contain the component 'functions' -- fit using the option model=TRUE")
  }

  ################################################################
  warning("Sorry, the profile likelihood is not implemented yet.")
  ###############################################################

  ## Control of the optim procedure
  if (is.null(control[["fnscale",exact=TRUE]])) { control$fnscale <- -1 }
  if (is.null(control[["maxit",exact=TRUE]])) { control$maxit <- 100 }
  if (is.null(control[["trace",exact=TRUE]])) { control$trace <- 2 }


  ## Estimated normalized likelihood function
  ltildeestim <- function(thetai,i) {
    theta <- theta.ml
    theta[i] <- thetai
    fitted$functions$ll(theta) - loglik.theta.ml
  }

  ## Profile normalized likelihood function
  ltildeprofile <- function(thetai,i)
  {
    #cat("Investigating theta[",i,"] = ",thetai,"\n")

    emptyTheta <- rep(0, length(theta.ml))

    # Likelihood l(theta_{-i}) = l(theta_i, theta_i)
    ltildethetaminusi <- function(thetaminusi) {
      theta <- emptyTheta
      theta[-i] <- thetaminusi
      theta[i] <- thetai
      #cat("Investigating theta = ",theta,"\n")
      res <- fitted$functions$ll(theta) - loglik.theta.ml
      #cat("Current ltildethetaminusi value: ",res,"\n")
      return(res)
    }
    # Score function of all params except thetaminusi
    stildethetaminusi <- function(thetaminusi) {
      theta <- emptyTheta
      theta[-i] <- thetaminusi
      theta[i] <- thetai
      res <- fitted$functions$sc(theta)[-i]
      #cat("Current stildethetaminusi value: ",res,"\n")
      return(res)
    }

    # Call optim -- currently not adapted to arguments of control arguments
    # used in the fit
    resOthers <- tryCatch(
            optim(par=theta.ml[-i], fn = ltildethetaminusi, gr = stildethetaminusi,
                  method = "BFGS", control = control),
            warning = function(w) print(w), error = function(e) list(value=NA))
    resOthers$value
  }



  ## Initialize
  theta.ml <- coef(fitted)
  loglik.theta.ml <- c(logLik(fitted))
  se <- sqrt(diag(vcov(fitted)))
  resProfile <- list()


  ## Perform profile computations for all requested parameters
  cat("Evaluating the profile logliks on a grid...\n")
  for (i in 1:length(profile))
    {
    cat("i= ",i,"/",length(profile),"\n")
    #Index of the parameter in the theta vector
    idx <- profile[[i]][1]
    #If no borders are given use those from wald intervals (unconstrained)
    if (is.na(profile[[i]][2])) profile[[i]][2] <- theta.ml[idx] - 3*se[idx]
    if (is.na(profile[[i]][3])) profile[[i]][3] <- theta.ml[idx] + 3*se[idx]
    #Evaluate profile loglik on a grid (if requested)
    if (profile[[i]][4] > 0) {
      thetai.grid <- seq(profile[[i]][2],profile[[i]][3],length=profile[[i]][4])
      resProfile[[i]] <- matrix(NA, nrow = length(thetai.grid), ncol = 4L,
        dimnames = list(NULL, c("grid","profile","estimated","wald")))

      #Loop over all gridpoints
      for (j in 1:length(thetai.grid)) {
        cat("\tj= ",j,"/",length(thetai.grid),"\n")
        resProfile[[i]][j,] <- c(thetai.grid[j],
           #Do we need to compute ltildeprofile (can be quite time consuming)
           ifelse(do.ltildeprofile, ltildeprofile(thetai.grid[j],idx), NA),
           ltildeestim(thetai.grid[j],idx),
           - 1/2*(1/se[idx]^2)*(thetai.grid[j] - theta.ml[idx])^2)
      }
    }
  }
  names(resProfile) <- names(theta.ml)[sapply(profile, function(x) x[1L])]

  ###############################
  ## Profile likelihood intervals
  ###############################
  # Not done, yet
  ciProfile <- NULL

  ####Done, return
  return(list(lp=resProfile, ci.hl=ciProfile, profileObj=profile))
}



### update-method for the twinstim-class
## stats::update.default would also work but is not adapted to the specific
## structure of twinstim: optim.args (use modifyList), two formulae, model, ...
## However, this specific method is inspired by and copies small parts of the
## update.default method from the stats package developed by The R Core Team

update.twinstim <- function (object, endemic, epidemic,
                             control.siaf, optim.args, model,
                             ..., use.estimates = TRUE, evaluate = TRUE)
{
    call <- object$call
    thiscall <- match.call(expand.dots=FALSE)
    extras <- thiscall$...
    
    if (!missing(model)) {
        call$model <- model
        ## Special case: update model component ONLY
        if (evaluate &&
            all(names(thiscall)[-1] %in% c("object", "model", "evaluate"))) {
            return(.update.twinstim.model(object, model))
        }
    }

    ## Why we no longer use call$endemic but update object$formula$endemic:
    ## call$endemic would be an unevaluated expression eventually receiving the
    ## parent.frame() as environment, cp.: 
    ##(function(e) {ecall <- match.call()$e; eval(call("environment", ecall))})(~1+start)
    ## This could cause large files if the fitted model is saved.
    ## Furthermore, call$endemic could refer to some object containing
    ## the formula, which is no longer visible.    
    call$endemic <- if (missing(endemic)) object$formula$endemic else
        update.formula(object$formula$endemic, endemic)
    call$epidemic <- if (missing(epidemic)) object$formula$epidemic else
        update.formula(object$formula$epidemic, epidemic)
    ## Note: update.formula uses terms.formula(...,simplify=TRUE), but
    ##       the principle order of terms is retained. Offsets will be moved to 
    ##       the end and a missing intercept will be denoted by a final -1.
    
    if (!missing(control.siaf))
        call$control.siaf <- if (is.null(object$control.siaf)) # if constantsiaf
            control.siaf else modifyList(object$control.siaf, control.siaf)
    
    call["optim.args"] <- if (missing(optim.args)) object["optim.args"] else {
        list( # use list() to enable optim.args=NULL
             if (is.list(optim.args)) {
                 modifyList(object$optim.args, optim.args)
             } else optim.args           # = NULL
             )
    }
    ## Set initial values (will be appropriately subsetted and/or extended with
    ## zeroes inside twinstim())
    call$start <- if (missing(optim.args) ||
                      (!is.null(optim.args) && !"par" %in% names(optim.args))) {
        ## old optim.args$par probably doesn't match updated model,
        ## thus we set it as "start"-argument
        call$optim.args$par <- NULL
        if (use.estimates) coef(object) else object$optim.args$par
    } else NULL
    if ("start" %in% names(extras)) {
        newstart <- eval.parent(extras$start)
        call$start[names(newstart)] <- newstart
        extras$start <- NULL
    }
    ## CAVE: the remainder is copied from stats::update.default (as at R-2.15.0)
    if(length(extras)) {
	existing <- !is.na(match(names(extras), names(call)))
	## do these individually to allow NULL to remove entries.
	for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	if(any(!existing)) {
	    call <- c(as.list(call), extras[!existing])
	    call <- as.call(call)
	}
    }
    if(evaluate) eval(call, parent.frame())
    else call
}

.update.twinstim.model <- function (object, model)
{
    call <- object$call
    call$model <- model
    if (model) { # add model environment
        call$start <- coef(object)
        call$optim.args$fixed <- TRUE
        call$cumCIF <- FALSE
        cat("Setting up the model environment ...\n")
        capture.output(suppressMessages(
                       objectWithModel <- eval(call, parent.frame())
                       ))
        cat("Done.\n")
        object <- append(object, objectWithModel["functions"],
                         after=match("optim.args", names(object)))
        class(object) <- class(objectWithModel) # dropped by append()
        environment(object) <- environment(objectWithModel)
    } else { # remove model environment
        object$functions <- NULL
        environment(object) <- NULL
    }
    object$call$model <- model
    object
}

## a terms-method is required for stepComponent()
terms.twinstim <- function (x, component=c("endemic", "epidemic"), ...)
{
    component <- match.arg(component)
    terms.formula(x$formula[[component]], keep.order=TRUE)
}

