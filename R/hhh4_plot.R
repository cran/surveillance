################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Plot-method(s) for fitted hhh4() models
###
### Copyright (C) 2010-2014 Michaela Paul and Sebastian Meyer
### $Revision: 833 $
### $Date: 2014-03-12 00:46:56 +0100 (Wed, 12 Mar 2014) $
################################################################################


plot.hhh4 <- function (x, type=c("fitted", "season", "maxEV", "ri"), ...)
{
    cl <- sys.call()
    cl$type <- NULL
    cl[[1L]] <- as.name(paste("plotHHH4", match.arg(type), sep="_"))
    eval(cl, envir=parent.frame())
}


###
### Time series of fitted component means and observed counts for selected units
###

plotHHH4_fitted <- function (x, units = 1, names = NULL,
                             col = c("orange","blue","grey85","black"),
                             pch = 19, pt.cex = 0.6,
                             par.settings = list(),
                             legend = TRUE, legend.args = list(),
                             legend.observed = TRUE, ...)
{
    if (!is.null(names)) stopifnot(length(units) == length(names))
    if (length(col) == 3) col <- c(col, "black") else if (length(col) != 4)
        stop("'col' must be of length 3 or 4") # for backwards compatibility

    if (is.list(par.settings)) {
        par.defaults <- list(mfrow = sort(n2mfrow(length(units))),
                             mar = c(4,4,2,0.5)+.1, las = 1)
        par.settings <- modifyList(par.defaults, par.settings)
        opar <- do.call("par", par.settings)
        on.exit(par(opar))
    }

    ## legend options
    if (is.logical(legend)) legend <- which(legend)
    if (!is.list(legend.args)) {
        if (length(legend) > 0)
            warning("ignored 'legend' since 'legend.args' is not a list")
        legend <- integer(0L)
    }
    if (length(legend) > 0) {
        legendidx <- 1L + c(if (legend.observed && !is.na(pch)) 0L,
                            which(c("ne","ar","end") %in% componentsHHH4(x)))
        default.args <- list(
            x="topright", col=c(col[4],col[1:3])[legendidx], lwd=6,
            lty=c(NA,1,1,1)[legendidx], pch=c(pch,NA,NA,NA)[legendidx],
            pt.cex=pt.cex, pt.lwd=1, bty="n", inset=0.02,
            legend=c("observed","spatiotemporal","autoregressive","endemic")[legendidx]
            )
        legend.args <- modifyList(default.args, legend.args)
    }

    ## plot fitted values region by region
    meanHHHunits <- vector(mode="list", length=length(units))
    for(i in seq_along(units)) {
        meanHHHunits[[i]] <- plotHHH4_fitted1(x, units[i], main=names[i],
                                              col=col, pch=pch, pt.cex=pt.cex,
                                              ...)
        if (i %in% legend) do.call("legend", args=legend.args)
    }
    invisible(meanHHHunits)
}


### plot estimated component means for a single region

plotHHH4_fitted1 <- function(x, unit=1, main=NULL,
                             col=c("grey30","grey60","grey85","grey0"),
                             pch=19, pt.cex=0.6, border=col,
                             start=x$stsObj@start, end=NULL,
                             xlim=NULL, ylim=NULL, xlab="", ylab="No. infected",
                             hide0s=FALSE, meanHHH=NULL)
{
    stsObj <- x$stsObj
    if (is.character(unit) &&
        is.na(unit <- match(.unit <- unit, colnames(stsObj))))
        stop("region '", .unit, "' does not exist")
    if (is.null(main)) main <- colnames(stsObj)[unit]

    ## get observed counts
    obs <- observed(stsObj)[,unit]

    ## time range for plotting
    start0 <- yearepoch2point(stsObj@start, stsObj@freq, toleft=TRUE)
    start <- yearepoch2point(start, stsObj@freq)
    tp <- start0 + seq_along(obs)/stsObj@freq # all observation time points
    if (start < start0 || start > tp[length(tp)])
        stop("'start' is not within the time range of 'x$stsObj'")
    end <- if(is.null(end)) tp[length(tp)] else yearepoch2point(end,stsObj@freq)
    stopifnot(start < end)
    tpInRange <- which(tp >= start & tp <= end)            # plot only those
    tpSubset <- tp[intersect(x$control$subset, tpInRange)] # fitted time points

    ## get fitted component means
    if (is.null(meanHHH))
        meanHHH <- meanHHH(coef(x,reparamPsi=FALSE), terms(x))
    meanHHHunit <- sapply(meanHHH, "[", i=TRUE, j=unit)
    stopifnot(is.matrix(meanHHHunit), !is.null(colnames(meanHHHunit)),
              nrow(meanHHHunit) == length(x$control$subset))
    meanHHHunit <- meanHHHunit[x$control$subset %in% tpInRange,,drop=FALSE]
    
    ## establish basic plot window
    if (is.null(ylim)) ylim <- c(0, max(obs[tpInRange],na.rm=TRUE))
    plot(c(start,end), ylim, xlim=xlim, xlab=xlab, ylab=ylab, type="n")
    title(main=main, line=0.5)

    ## draw polygons
    xpoly <- c(tpSubset[1], tpSubset, tail(tpSubset,1))
    polygon(xpoly, c(0,meanHHHunit[,"mean"],0),
            col=col[1], border=border[1])
    if (x$control$ar$inModel)
        polygon(xpoly, c(0,rowSums(meanHHHunit[,c("endemic","epi.own")]),0),
                col=col[2], border=border[2])
    if (x$control$end$inModel)
        polygon(xpoly, c(0,meanHHHunit[,"endemic"],0),
                col=col[3], border=border[3])

    ## add observed counts within [start;end]
    ptidx <- if (hide0s) intersect(tpInRange, which(obs > 0)) else tpInRange
    points(tp[ptidx], obs[ptidx], col=col[4], pch=pch, cex=pt.cex)

    ## invisibly return the fitted component means for the selected region
    invisible(meanHHHunit)
}


###
### Map of estimated random intercepts of a specific component
###

plotHHH4_ri <- function (x, component, sp.layout = NULL,
                         gpar.missing = list(col="darkgrey", lty=2, lwd=2),
                         ...)
{
    ranefmatrix <- ranef(x, tomatrix=TRUE)
    if (is.null(ranefmatrix)) stop("model has no random effects")
    stopifnot(length(component) == 1L)
    if (is.na(comp <- pmatch(component, colnames(ranefmatrix))))
        stop("'component' must (partially) match one of ",
             paste(dQuote(colnames(ranefmatrix)), collapse=", "))
    
    map <- x$stsObj@map
    if (is.null(map)) stop("'x$stsObj' has no map")
    map$ranef <- ranefmatrix[,comp][row.names(map)]
    
    if (is.list(gpar.missing) && any(is.na(map$ranef)))
        sp.layout <- c(sp.layout, 
                       c(list("sp.polygons", map[is.na(map$ranef),]),
                         gpar.missing))
    spplot(map[!is.na(map$ranef),], zcol = "ranef",
           sp.layout = sp.layout, ...)
}


###
### Plot the course of the dominant eigenvalue of one or several hhh4-fits
###

plotHHH4_maxEV <- function (...,
                            matplot.args = list(), refline.args = list(),
                            legend.args = list())
{
    objnams <- unlist(lapply(match.call(expand.dots=FALSE)$..., deparse))
    objects <- getHHH4list(..., .names = objnams)

    ## get time points
    epoch <- attr(objects, "epoch")
    start <- attr(objects, "start")
    freq <- attr(objects, "freq")
    start0 <- yearepoch2point(start, freq, toleft=TRUE)
    tp <- start0 + seq_along(epoch) / freq

    ## compute course of dominant eigenvalue for all models
    maxEV <- sapply(objects, getMaxEV, simplify=TRUE, USE.NAMES=TRUE)

    ## line style
    matplot.args <- modifyList(
        list(type="l", col=c(1,2,6,3), lty=c(1,3,2,4), lwd=1.7, cex=1, pch=NULL,
             xlab="", ylab="dominant eigenvalue", ylim=c(0,max(2,maxEV))),
        matplot.args)

    ## main plot
    do.call("matplot", c(list(x=tp, y=maxEV), matplot.args))

    ## add reference line
    if (is.list(refline.args))
        do.call("abline", modifyList(list(h=1, lty=3, col="grey"),
                                     refline.args))

    ## add legend
    if (missing(legend.args) && length(objects) == 1)
        legend.args <- NULL             # omit legend
    if (is.list(legend.args)) {
        legend.args <- modifyList(
            c(list(x="topright", inset=0.02, legend=names(objects), bty="n"),
              matplot.args[c("col", "lwd", "lty", "pch")],
              with(matplot.args, list(pt.cex=cex, text.col=col))),
            legend.args)
        do.call("legend", legend.args)
    }

    ## done
    invisible(maxEV)
}

getMaxEV <- function (x)
{
    Lambda <- createLambda(x)
    if (identical(type <- attr(Lambda, "type"), "zero")) {
        rep.int(0, nrow(x$stsObj))
    } else {
        sapply(seq_len(nrow(x$stsObj)),
               if (identical(type, "diagonal")) {
                   function (t) max(Lambda(t))  # or max(diag(Lambda(t)))
               } else { # type is NULL
                   function (t)
                       eigen(Lambda(t), symmetric=FALSE,
                             only.values=TRUE)$values[1L]
               }, simplify=TRUE, USE.NAMES=FALSE)
    }
}

## generate a function that computes the Lambda_t matrix
createLambda <- function (object)
{
    nTime <- nrow(object$stsObj)
    nUnit <- object$nUnit
    if (identical(componentsHHH4(object), "end")) { # no epidemic components
        zeromat <- matrix(0, nUnit, nUnit)
        Lambda <- function (t) zeromat
        attr(Lambda, "type") <- "zero"
        return(Lambda)
    }
    
    meanHHH <- meanHHH(object$coefficients, terms(object),
                       subset=seq_len(nTime))
    
    W <- getNEweights(object)
    Wt <- if (is.null(W)) {
        NULL
    } else if (is.matrix(W)) {
        function (t) W
    } else {
        function (t) W[,,t]
    }

    type <- NULL
    Lambda <- if (is.null(Wt)) {        # no neighbourhood component
        type <- "diagonal"
        function (t) {
            stopifnot(isScalar(t) && t > 0 && t <= nTime)
            diag(meanHHH$ar.exppred[t,], nUnit, nUnit)
        }
    } else {
        function (t) {
            stopifnot(isScalar(t) && t > 0 && t <= nTime)
            Lambda <- meanHHH$ne.exppred[t,] * Wt(t)
            diag(Lambda) <- meanHHH$ar.exppred[t,]
            Lambda
        }
    }
    attr(Lambda, "type") <- type
    Lambda
}


###
### Plot estimated seasonality (sine-cosine terms) of one or several hhh4-fits
###

plotHHH4_season <- function (...,
                             components = c("ar", "ne", "end", "maxEV"),
                             xlim = NULL, ylim = NULL,
                             xlab = NULL, ylab = NULL, main = NULL,
                             par.settings = list(), matplot.args = list(),
                             legend = NULL, legend.args = list(),
                             refline.args = list(), unit = 1)
{
    objnams <- unlist(lapply(match.call(expand.dots=FALSE)$..., deparse))
    objects <- getHHH4list(..., .names = objnams)
    freq <- attr(objects, "freq")
    components <- match.arg(components, several.ok=TRUE)

    ## x-axis
    if (is.null(xlim))
        xlim <- c(1,freq)
    if (is.null(xlab))
        xlab <- if (freq==52) "week" else if (freq==12) "month" else "time"

    ## auxiliary function for an argument list "x" with named "defaults" list
    withDefaults <- function(x, defaults)
    {
        if (is.null(x)) defaults else if (is.list(x)) {
            if (is.null(names(x))) {    # x must be complete
                stopifnot(length(x) == length(defaults))
                setNames(x, names(defaults))
            } else modifyList(defaults, x) # x might be a subset of parameters
        } else if (is.atomic(x)) {
            setNames(rep(list(x), length(defaults)), names(defaults))
        } else stop("'", deparse(substitute(x)), "' is not suitably specified")
    }
    
    ## component-specific arguments 
    ylim <- withDefaults(ylim,
                         list(ar=NULL, ne=NULL, end=NULL, maxEV=NULL))
    ylab <- withDefaults(ylab,
                         list(ar=expression(hat(lambda)),
                              ne=expression(hat(phi)),
                              end=expression(hat(nu)),
                              maxEV="dominant eigenvalue"))
    main <- withDefaults(main,
                         list(ar="autoregressive component",
                              ne="spatiotemporal component",    
                              end="endemic component",
                              maxEV="dominant eigenvalue"))
    anyMain <- any(unlist(lapply(main, nchar),
                          recursive=FALSE, use.names=FALSE) > 0)

    ## basic graphical settings
    if (is.list(par.settings)) {
        par.defaults <- list(mfrow=sort(n2mfrow(length(components))),
                             mar=c(4,5,if(anyMain) 2 else 1,1)+.1, las=1)
        par.settings <- modifyList(par.defaults, par.settings)
        opar <- do.call("par", par.settings)
        on.exit(par(opar))
    }

    ## line style
    matplot.args <- modifyList(list(type="l", col=c(1,2,6,3), lty=c(1,3,2,4),
                                    lwd=1.7, cex=1, pch=NULL),
                               matplot.args)

    ## legend options
    if (is.null(legend)) legend <- length(objects) > 1
    if (is.logical(legend)) legend <- which(legend)
    if (!is.list(legend.args)) {
        if (length(legend) > 0)
            warning("ignored 'legend' since 'legend.args' is not a list")
        legend <- integer(0L)
    }
    if (length(legend) > 0) {
        default.args <- c(
            list(x="topright", inset=0.02, legend=names(objects), bty="n"),
            matplot.args[c("col", "lwd", "lty", "pch")],
            with(matplot.args, list(pt.cex=cex, text.col=col))
            )
        legend.args <- modifyList(default.args, legend.args)
    }

    ## plot seasonality in individual model components
    seasons <- list()
    for(comp in setdiff(components, "maxEV")){
        s2 <- lapply(objects, getSeason, component = comp, unit = unit)
        seasons[[comp]] <- exp(sapply(s2, function(intseas) do.call("+", intseas)))
        do.call("matplot",              # x defaults to 1:freq
                c(list(seasons[[comp]], xlim=xlim,
                       ylim=if (is.null(ylim[[comp]]))
                       c(0,max(1,seasons[[comp]])) else ylim[[comp]],
                       xlab=xlab, ylab=ylab[[comp]],
                       main=main[[comp]]), matplot.args))
        if (match(comp, components) %in% legend)
            do.call("legend", legend.args)
    }

    ## plot seasonality of dominant eigenvalue
    if ("maxEV" %in% components) {
        seasons[["maxEV"]] <- sapply(objects, function(obj)
                                     getMaxEV_season(obj, unit=unit)$maxEV.season)
        do.call("matplot",
                c(list(seasons[["maxEV"]], xlim=xlim,
                       ylim=if (is.null(ylim[["maxEV"]]))
                       c(0,max(2,seasons[["maxEV"]])) else ylim[["maxEV"]],
                       xlab=xlab, ylab=ylab[["maxEV"]],
                       main=main[["maxEV"]]), matplot.args))
        if (is.list(refline.args))
            do.call("abline", modifyList(list(h=1, lty=3, col="grey"),
                                         refline.args))
        if (4 %in% legend) do.call("legend", legend.args)
    }

    ## invisibly return the data that has been plotted
    invisible(seasons)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get estimated seasonal pattern in the different components
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getSeason <- function(x, component = c("ar", "ne", "end"), unit = 1)
{
    stopifnot(inherits(x, c("hhh4","ah4")))
    component <- match.arg(component)
    startseason <- getSeasonStart(x)
    freq <- x$stsObj@freq
    if (is.character(unit)) unit <- match(unit, colnames(x$stsObj))
    
    ## return -Inf is component is not in the model (-> exp(-Inf) = 0)
    if (!component %in% componentsHHH4(x))
        return(list(intercept=-Inf, season=rep.int(-Inf, freq)))

    ## get the intercept
    est <- fixef(x, reparamPsi=FALSE)
    intercept <- est[grep(paste0("^", component, "\\.(1|ri)"), names(est))]
    if (length(intercept) == 0) {
        intercept <- 0 # no intercept (not standard)
    } else if (length(intercept) > 1) { # unit-specific intercepts
        if (length(intercept) != ncol(x$stsObj))
            stop(component,"-component has incomplete unit-specific intercepts")
        intercept <- intercept[unit]
        if (is.na(intercept)) stop("the specified 'unit' does not exist")
    }

    ## get seasonality terms (relying on sin(2*pi*t/52)-kind coefficient names)
    coefSinCos <- est[grep(paste0("^",component, "\\.(sin|cos)\\("), names(est))]
    if (unitspecific <- length(grep(").", names(coefSinCos), fixed=TRUE))) {
        if (unitspecific < length(coefSinCos))
            stop("cannot handle partially unit-specific seasonality")
        coefSinCos <- coefSinCos[grep(paste0(").",colnames(x$stsObj)[unit]),
                                      names(coefSinCos), fixed=TRUE)]
    }
    if (length(coefSinCos)==0)
        return(list(intercept=intercept, season=rep.int(0,freq)))
    fSinCos <- reformulate(
        sub(paste0("^",component,"\\."), "", names(coefSinCos)),
        intercept=FALSE)
    mmSinCos <- model.matrix(fSinCos,
                             data=data.frame(t=startseason-1 + seq_len(freq)))

    ## Done
    list(intercept=intercept, season=mmSinCos %*% coefSinCos)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compute dominant eigenvalue of Lambda_t
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getMaxEV_season <- function (x, unit = 1)
{
    stopifnot(inherits(x, c("hhh4","ah4")))
    nUnits <- x$nUnit
    freq <- x$stsObj@freq
    components <- componentsHHH4(x)
    
    ## global intercepts and seasonality
    s2.lambda <- getSeason(x, "ar", unit = unit)
    s2.phi <- getSeason(x, "ne", unit = unit)
    
    ## unit-specific intercepts
    ris <- ranef(x, tomatrix=TRUE)
    ri.lambda <- ris[,pmatch("ar.ri", colnames(ris), nomatch=0L),drop=TRUE]
    if (length(ri.lambda) == 0L) ri.lambda <- rep.int(0, nUnits)
    ri.phi <- ris[,pmatch("ne.ri", colnames(ris), nomatch=0L),drop=TRUE]
    if (length(ri.phi) == 0L) ri.phi <- rep.int(0, nUnits)

    ## get neighbourhood weights as a function of time
    W <- getNEweights(x)  # NULL, matrix or 3-dim array
    if (!is.null(W) && !is.matrix(W))
        stop("neighbourhood weights are time-varying; ",
            # and thus probably changing within or across seasons
             "use getMaxEV() or plotHHH4_maxEV()")

    ## create the Lambda_t matrix
    createLambda <- function (t)
    {
        Lambda <- if ("ne" %in% components) {
            exp(s2.phi$intercept + ri.phi + if(t==0) 0 else s2.phi$season[t]) * W
        } else matrix(0, nUnits, nUnits)
        diag(Lambda) <- if ("ar" %in% components) {
            exp(s2.lambda$intercept + ri.lambda +
                if(t==0) 0 else s2.lambda$season[t])
        } else 0
        Lambda
    }

    ## calculate the maximum eigenvalue of Lambda_t
    maxEV <- function(t)
    {
        ev <- eigen(createLambda(t), only.values=TRUE)$values
        if (is.complex(ev)) {
            warning("Lambda_",t," has complex eigenvalues")
            ev <- abs(ev)   # if complex, use abs EVs
        }
        max(ev)
    }

    ## do this for t in 0:freq
    maxEV.const <- maxEV(0)
    maxEV.season <- if (all(c(s2.phi$season, s2.lambda$season) %in% c(-Inf, 0)))
        rep.int(maxEV.const, freq) else sapply(seq_len(freq), maxEV)

    ## Done
    list(maxEV.season = maxEV.season,
         maxEV.const  = maxEV.const,
         Lambda.const = createLambda(0))
}


## Determine the time point t of the start of a season in a hhh4() fit.
## If \code{object$stsObj@start[2] == 1}, it simply equals
## \code{object$control$data$t[1]}. Otherwise, the \code{stsObj} time series
## starts within a year (at sample \code{s}, say) and the beginning of
## the next season is
## \code{object$control$data$t[1] + object$stsObj@freq - s + 1}.
getSeasonStart <- function (object)
{
    if ((startsample <- object$stsObj@start[2]) == 1) {
        object$control$data$t[1]
    } else {
        object$control$data$t[1] + object$stsObj@freq-startsample + 1
    }
}



### auxiliary functions

yearepoch2point <- function (yearepoch, frequency, toleft=FALSE)
    yearepoch[1L] + (yearepoch[2L] - toleft) / frequency


getHHH4list <- function (..., .names = NA_character_)
{
    objects <- list(...)
    if (length(objects) == 1L && is.list(objects[[1L]]) &&
        inherits(objects[[1L]][[1L]], c("hhh4","ah4"))) {
        ## ... is a single list of fits
        objects <- objects[[1L]]
        if (is.null(names(objects))) names(objects) <- seq_along(objects)
    } else {
        names(objects) <- if (is.null(names(objects))) .names else {
            ifelse(nzchar(names(objects)), names(objects), .names)
        }
    }
    if (!all(sapply(objects, inherits, what=c("hhh4","ah4"))))
        stop("'...' must consist of hhh4()-fits only")

    ## check common epoch, start and frequency and append them as attributes
    epoch <- unique(t(sapply(objects, function(x) x$stsObj@epoch)))
    if (nrow(epoch) > 1)
        stop("supplied hhh4-models obey different 'epoch's")
    attr(objects, "epoch") <- drop(epoch)
    start <- unique(t(sapply(objects, function(x) x$stsObj@start)))
    if (nrow(start) > 1)
        stop("supplied hhh4-models obey different start times")
    attr(objects, "start") <- drop(start)
    freq <- unique(sapply(objects, function(x) x$stsObj@freq))
    if (length(freq)>1)
        stop("supplied hhh4-models obey different frequencies")
    attr(objects, "freq") <- freq
    
    ## done
    return(objects)
}
