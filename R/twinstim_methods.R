################################################################################
### Methods for objects of class "twinstim", specifically:
### coef (coefficients), vcov, logLik, print, summary, print.summary, plot
### Author: Sebastian Meyer
################################################################################

### don't need a specific coef-method (identical to stats:::coef.default)
## coef.twinstim <- function (object, ...)
## {
##     object$coefficients
## }

# asymptotic variance-covariance matrix (inverse of fisher information matrix)
vcov.twinstim <- function (object, ...)
{
    solve(object$fisherinfo)  # inverse of estimated expected fisher information
}

logLik.twinstim <- function (object, ...)
{
    r <- object$loglik
    attr(r, "df") <- length(coef(object))
    class(r) <- "logLik"
    r
}

print.twinstim <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")
    print.default(x$call)
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits=digits), print.gap = 2, quote = FALSE)
    cat("\nLog-likelihood: ", format(logLik(x), digits=digits), "\n", sep = "")
    if (!x$converged) {
        cat("\nWARNING: OPTIMIZATION ROUTINE DID NOT CONVERGE!\n")
    }
    cat("\n")
    invisible(x)
}

summary.twinstim <- function (object, test.iaf = FALSE,
    correlation = FALSE, symbolic.cor = FALSE, ...)
{
    ans <- object[c("call", "converged", "counts")]
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
    if (q > 0L) {
        cat("\nCoefficients of the epidemic component:\n")
        printCoefmat(rbind(x$coefficients.gamma, x$coefficients.iaf), digits = digits,
            signif.stars = signif.stars, ...)
#         if (niafpars > 0L) {
#             #cat("Coefficients of interaction functions:\n")
#             printCoefmat(x$coefficients.iaf, digits = digits, signif.stars = signif.stars, ...)
#         }
    } else cat("\nNo epidemic component.\n")
#     nEvents <- x$nEvents
#     nh0 <- length(nEvents)
#     if (nh0 < 2L) {
#         cat("\nTotal number of infections: ", nEvents, "\n")
#     } else {
#         cat("\nBaseline intervals:\n")
#         intervals <- character(nh0)
#         for(i in seq_len(nh0)) {
#             intervals[i] <-
#             paste("(",
#                   paste(format(x$intervals[c(i,i+1L)],trim=TRUE), collapse=";"),
#                   "]", sep = "")
#         }
#         names(intervals) <- paste("logbaseline", seq_len(nh0), sep=".")
#         print.default(rbind("Time interval" = intervals,
#                             "Number of events" = nEvents),
#                       quote = FALSE, print.gap = 2)
#     }
#     cat("\n", attr(x$aic, "type"), ": ", format(x$aic, digits=max(4, digits+1)),
#         if (!attr(x$aic, "exact")) "\t(simulated penalty weights)" else "",
#         sep = "")
    cat("\nAIC: ", format(x$aic, digits=max(4, digits+1)))
    cat("\nLog-likelihood:", format(x$loglik, digits = digits))
    cat("\nNumber of log-likelihood evaluations:", x$counts[1])
    cat("\nNumber of score function evaluations:", x$counts[2])
    if (!is.null(x$runtime)) {
        cat("\nRuntime:", format(x$runtime, digits=max(4, digits+1)), "seconds")
    }
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
            correl <- format(round(correl, 2), nsmall = 2, digits = digits)
            correl[!lower.tri(correl)] <- ""
            print(correl[-1, -p, drop = FALSE], quote = FALSE)
        }
        }
    }
    if (!x$converged) {
        cat("\nWARNING: OPTIMIZATION ROUTINE DID NOT CONVERGE!\n")
    }
    cat("\n")
    invisible(x)
}



### 'cat's the summary in LaTeX code

toLatex.summary.twinstim <- function (object, digits = max(3, getOption("digits") - 3),
                                      align = "rrrrr", withAIC = TRUE, ...)
{
ret <- capture.output({
    cat("\\begin{tabular}{", align, "}\n\\hline\n", sep="")
    cat(" & Estimate & Std. Error & $z$ value & $\\P(|Z|>|z|)$ \\\\\n\\hline\n\\hline\n")

    tabh <- object$coefficients.beta
    tabe <- rbind(object$coefficients.gamma, object$coefficients.iaf)
    for (tabname in c("tabh", "tabe")) {
        tab <- get(tabname)
        if (nrow(tab) > 0L) {
            rownames(tab) <- gsub(" ", "", rownames(tab))
            tab_char <- capture.output(
                        printCoefmat(tab,digits=digits,signif.stars=FALSE,na.print="NA")
                        )[-1]
            tab_char <- sub("([<]?)[ ]?([0-9]+)e([+-][0-9]+)$",
                            "\\1\\2\\\\cdot{}10^{\\3}",
                            tab_char)
            con <- textConnection(tab_char)
            tab2 <- read.table(con, colClasses="character")
            close(con)
            rownames(tab2) <- paste0("\\texttt{",tab2[,1],"}")
            tab2 <- tab2[,-1]
            tab2[] <- lapply(tab2, function(x) {
                ifelse(is.na(x), "", paste0("$",x,"$")) # (test.iaf=FALSE)
            })
            print(xtable::xtable(tab2), only.contents=TRUE, hline.after=NULL,
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


### Plot temporal or spatial evolution of the intensity

intensityplot.twinstim <- function (x,
    which = c("epidemic proportion", "endemic proportion", "total intensity"),
    aggregate = c("time", "space"),
    types = 1:nrow(x$qmatrix), tiles, tiles.idcol = NULL,
    plot = TRUE, add = FALSE, tgrid = 101, rug.opts = list(),
    sgrid = 128, polygons.args = list(), points.args = list(cex=0.5),
    cex.fun = sqrt, ...)
{
    ## check arguments
    if (is.null(environment(x))) {
        stop("'", substitute(x), "' is missing the model environment -- re-fit with 'model=TRUE'")
    }
    which <- match.arg(which)
    aggregate <- match.arg(aggregate)
    stopifnot(is.vector(types, mode="numeric"), types %in% 1:nrow(x$qmatrix),
              !anyDuplicated(types))

    ## set up desired intensities
    cl <- match.call()
    cl <- cl[c(1L, match(names(formals(intensity.twinstim)), names(cl), 0L))]
    cl[[1]] <- as.name("intensity.twinstim")
    components <- eval(cl, envir = parent.frame())

    ## define function to plot
    FUN <- function (tmp) {}
    names(formals(FUN)) <- if (aggregate == "time") "times" else "coords"
    body1 <- if (aggregate == "time") expression(
        hGrid <- sapply(times, components$hFUN, USE.NAMES=FALSE),
        eGrid <- sapply(times, components$eFUN, USE.NAMES=FALSE)
        ) else expression(
            hGrid <- unname(components$hFUN(coords)), # works with whole coord matrix
            eGrid <- apply(coords, 1, components$eFUN)
            )
    body2 <- switch(which,
                    "epidemic proportion" = expression(eGrid / (hGrid + eGrid)),
                    "endemic proportion" = expression(hGrid / (hGrid + eGrid)),
                    "total intensity" = expression(hGrid + eGrid))
    body(FUN) <- as.call(c(as.name("{"), c(body1, body2)))
    
    if (!plot) return(FUN)

    ## plot the FUN
    modelenv <- environment(x)
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
            if (!require("maptools")) {
                stop("auto-generation of 'sgrid' requires package \"maptools\"")
            }
            sgrid <- maptools::Sobj_SpatialGrid(.tiles, n = sgrid)$SG
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
            nms.points <- names(points.args)
            if(! "pch" %in% nms.points) points.args$pch <- 1
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

intensity.twinstim <- function (x, aggregate = c("time", "space"),
    types = 1:nrow(x$qmatrix), tiles, tiles.idcol = NULL)
{
    ## check arguments
    if (is.null(environment(x))) {
        stop("'x' is missing the model environment -- re-fit with 'model=TRUE'")
    }
    aggregate <- match.arg(aggregate)
    stopifnot(is.vector(types, mode="numeric"), types %in% 1:nrow(x$qmatrix),
              !anyDuplicated(types))

    ## model environment
    modelenv <- environment(x)          # for the functions to be returned
    thisenv <- environment()
    parent.env(thisenv) <- modelenv     # objects of modelenv become visible
    ## CAVE: The R manual says:
    ## "The replacement function 'parent.env<-' is extremely dangerous [...].
    ##  It may be removed in the near future."
    qmatrix <- x$qmatrix                # not part of modelenv
    force(types)                        # evaluate types before rm(x)
    rm(x)                               # don't need this anymore
    
    ## endemic component on the spatial or temporal grid
    hInt <-
        if (hash) {
            eta <- drop(mmhGrid %*% beta)
            if (!is.null(offsetGrid)) eta <- offsetGrid + eta
            expeta <- exp(unname(eta))
            .beta0 <- rep(if (nbeta0==0L) 0 else beta0, length = nTypes)
            fact <- sum(exp(.beta0[types]))
            if (aggregate == "time") {      # int over W and types by BLOCK
                fact * c(tapply(expeta * ds, gridBlocks, sum, simplify = TRUE))
            } else {                        # int over T and types by tile
                fact * c(tapply(expeta * dt, gridTiles, sum, simplify = TRUE))
            }
        } else {
            ngrid <- if (aggregate == "time") {
                gridBlocks[length(gridBlocks)]
            } else nlevels(gridTiles)
            rep.int(0, ngrid)
        }

    ## endemic component as a function of time or location
    hIntFUN <- if (hash) {
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
                polygonidxOfPoints <- sp::over(points, .tiles)
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
    eInt <- if (hase) {
        qSum_types <- rowSums(qmatrix[,types,drop=FALSE])[eventTypes]
        fact <- qSum_types * gammapred
        if (aggregate == "time") {      # as a function of time (int over W & types)
            factS <- fact * siafInt
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
        } else {                        # as a function of location (int over time and types)
            factT <- fact * tiafInt
            function (xy) {
                stopifnot(is.vector(xy, mode="numeric"), length(xy) == 2L)
                point <- matrix(xy, nrow=nrow(eventCoords), ncol=2L, byrow=TRUE)
                sdiff <- point - eventCoords
                proximity <- qSum_types > 0 & rowSums(sdiff^2) <= eps.s^2
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
    list(hGrid = hInt, hFUN = hIntFUN, eFUN = eInt)
}


### Plot fitted tiaf or siaf(cbind(0, r)), r=distance

iafplot <- function (object, which = c("siaf", "tiaf"),
    types = 1:nrow(object$qmatrix),
    conf.type = if (length(pars) > 1) "bootstrap" else "parbounds",
    conf.level = 0.95, conf.B = 999,
    ngrid = 101, col.estimate = rainbow(length(types)), col.conf = col.estimate,
    alpha.B = 0.15, lwd = c(3,1), lty = c(1,2), xlim = c(0,eps), ylim = c(0,1),
    add = FALSE, xlab = NULL, ylab = NULL, ...)
{
    which <- match.arg(which)
    eps <- if (which == "siaf") {
        sqrt(sum((object$bbox[,"max"] - object$bbox[,"min"])^2))
    } else {
        diff(object$timeRange)
    }
    FUN <- object$formula[[which]][[if (which=="siaf") "f" else "g"]]
    coefs <- coef(object)
    idxpars <- grep(which,names(coefs))
    pars <- coefs[idxpars]
    force(conf.type)
    if (length(pars) == 0 || any(is.na(conf.type)) || is.null(conf.type)) {
        conf.type <- "none"
    }
    conf.type <- match.arg(conf.type, choices = c("parbounds", "bootstrap", "none"))
    
    if (!add) {
        if (is.null(xlab)) xlab <- if (which == "siaf") {
            expression("Distance " * group("||",bold(s)-bold(s)[j],"||") * " from host")
        } else {
            expression("Distance " * t-t[j] * " from host")
        }
        if (is.null(ylab)) ylab <- if (which == "siaf") {
            expression(f(group("||",bold(s)-bold(s)[j],"||")))
        } else {
            expression(g(t-t[j]))
        }
        plot(xlim, ylim, type="n", xlab = xlab, ylab = ylab, ...)
    }
    xgrid <- seq(0, xlim[2], length.out=ngrid)
    
    for (i in seq_along(types)) {
        ## select parameters on which to evaluate iaf
        parSample <- switch(conf.type, parbounds = {
            cis <- confint(object, idxpars, level=conf.level)
            ## all combinations of parameter bounds
            do.call("expand.grid", as.data.frame(t(cis)))
        }, bootstrap = {
            ## bootstrapping parameter values
            if (length(pars) == 1) {
                as.matrix(c(pars, rnorm(conf.B, mean=pars,
                                sd=sqrt(vcov(object)[idxpars,idxpars]))))
            } else {
                rbind(pars, mvtnorm::rmvnorm(conf.B, mean=pars,
                                 sigma=vcov(object)[idxpars,idxpars,drop=FALSE]))
            }
        })
        
        ## add confidence limits
        if (!is.null(parSample)) {
            fvalsSample <- apply(parSample, 1, function(pars)
                                 FUN(if(which=="siaf") cbind(xgrid,0) else xgrid,
                                     pars, types[i]))
            lowerupper <- if (conf.type == "parbounds") {
                t(apply(fvalsSample, 1, range))
            } else { # bootstrapped parameter values
                if (is.na(conf.level)) {
                    stopifnot(alpha.B >= 0, alpha.B <= 1)
                    .col <- col2rgb(col.conf[i], alpha=TRUE)[,1]
                    .col["alpha"] <- round(alpha.B*.col["alpha"])
                    .col <- do.call("rgb", args=c(as.list(.col), maxColorValue = 255))
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
                     " -- re-fit with 'model=TRUE'")
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
        if (attr(form$siaf, "constant")) {
            iRareas <- sapply(influenceRegion, spatstat::area.owin)
            ## will be used by .siafInt()
        } else if (! (is.null(siaf$Fcircle) ||
               (is.null(siaf$effRange) && noCircularIR))) {
            if (is.null(bdist)) {
                stop("missing \".bdist\" component in 'newevents'")
            }
        }
        nCub <- object$nCub
        .siafInt <- .siafIntFUN(siaf, nCub.adaptive=object$nCub.adaptive, noCircularIR=noCircularIR)
        siafInt <- .siafInt(siafpars)
        
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
            if (is.null(siaf$Fcircle)) { # implies that nCub was non-adaptive
                polyCub.midpoint(discpoly(c(0,0), eps.s, class="owin"),
                                 form$siaf$f, siafpars, type, eps=object$nCub)
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
## However, this specific method is inspired by and copies parts of the
## update.default method from the stats package developed by The R Core Team

update.twinstim <- function (object, endemic, epidemic, optim.args, model,
                             ..., evaluate = TRUE)
{
    call <- object$call
    thiscall <- match.call(expand.dots=FALSE)

    if (!missing(model)) {
        call$model <- model
        ## Special case: update model component only
        if (evaluate &&
            all(names(thiscall)[-1] %in% c("object", "model", "evaluate"))) {
            if (model) {                # add model environment
                call$optim.args$par <- coef(object)
                call$optim.args$fixed <- seq_along(coef(object))
                call$cumCIF <- FALSE
                cat("Setting up the model environment ...\n")
                capture.output(suppressMessages(
                    objectWithModel <- eval(call, parent.frame())
                ))
                cat("Done.\n")
                object$formula <- objectWithModel$formula
                object$functions <- objectWithModel$functions
                ## CAVE: order of components of the twinstim-object
                object <- object[c(1:match("formula", names(object)),
                                 match(c("functions", "call", "runtime"), names(object)))]
                environment(object) <- environment(objectWithModel)
                ## re-add class attribute (dropped on re-ordering)
                class(object) <- "twinstim"
            } else {                    # remove model environment
                environment(object$formula$epidemic) <-
                    environment(object$formula$endemic) <- .GlobalEnv
                object$functions <- NULL
                environment(object) <- NULL
            }
            object$call$model <- model
            return(object)
        }
    }
    if (!missing(endemic))
        call$endemic <- stats::update.formula(formula(object)$endemic, endemic)
    if (!missing(epidemic))
        call$epidemic <- stats::update.formula(formula(object)$epidemic, epidemic)
    if (!missing(optim.args)) {
        oldargs <- call$optim.args
        call$optim.args <-
            if (is.listcall(oldargs) && is.list(optim.args)) {
                modifyListcall(oldargs, thiscall$optim.args)
            } else thiscall$optim.args
    }
    extras <- thiscall$...
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
