################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Plots for an array "hhh4sims" of simulated counts from an "hhh4" model,
### or a list thereof as produced by different "hhh4" models (same period!)
###
### Copyright (C) 2013-2018 Sebastian Meyer
### $Revision: 2066 $
### $Date: 2018-01-22 11:46:23 +0100 (Mon, 22. Jan 2018) $
################################################################################

plot.hhh4sims <- function (x, ...) {
    ## use the object name of x
    x <- eval(substitute(as.hhh4simslist(x)), envir = parent.frame())
    plot.hhh4simslist(x, ...)
}


## class for a list of "hhh4sims" arrays from different models
## (over the same period with same initial values)
hhh4simslist <- function (x, initial, stsObserved)
{
    ## drop attributes from every single hhh4sims object
    for (i in seq_along(x))
        attr(x[[i]], "class") <- attr(x[[i]], "initial") <-
            attr(x[[i]], "stsObserved") <- NULL
    ## set as list attributes
    attr(x, "initial") <- initial
    attr(x, "stsObserved") <- stsObserved
    class(x) <- "hhh4simslist"
    x
}

## converter functions
as.hhh4simslist <- function (x, ...) UseMethod("as.hhh4simslist")

as.hhh4simslist.hhh4sims <- function (x, ...)
{
    ## we do not use x here, but construct a list() from the sys.call()
    ## such that as.hhh4simslist(name1 = model1, name2 = model2) works
    cl <- sys.call()
    cl[[1L]] <- as.name("list")
    xx <- eval(cl, envir = parent.frame())
    objnames <- as.character(cl)[-1L]
    if (is.null(names(xx))) {
        names(xx) <- objnames
    } else {
        names(xx)[names(xx) == ""] <- objnames[names(xx) == ""]
    }
    as.hhh4simslist.list(xx)
}

as.hhh4simslist.list <- function (x, ...)
{
    ## verify class
    lapply(X = x, FUN = function (Xi)
        if (!inherits(Xi, "hhh4sims"))
            stop(sQuote("x"), " is not a list of ", dQuote("hhh4sims")))
    hhh4simslist(x,
                 initial = attr(x[[1L]], "initial"),
                 stsObserved = attr(x[[1L]], "stsObserved"))
}

as.hhh4simslist.hhh4simslist <- function (x, ...) x


## 'x[i]': select models (elements of the list)
## 'x[i,j,]': subset simulations while keeping attributes in sync
"[.hhh4simslist" <- function (x, i, j, ..., drop = FALSE)
{
    ## case 1: select models
    if (nargs() == 2L) {
        ## select elements of the list
        xx <- NextMethod("[")
        ## restore class attributes
        xx <- hhh4simslist(xx,
                           initial = attr(x, "initial"),
                           stsObserved = attr(x, "stsObserved"))
        return(xx)
    }

    ## case 2: subset simulations, i.e., index individual arrays
    cl <- sys.call()
    cl[[1L]] <- as.name("[")
    cl[[2L]] <- quote(x)
    cl$drop <- drop
    subseti <- as.function(c(alist(x=), cl), envir = parent.frame())
    x[] <- lapply(X = unclass(x), subseti)  # unclass to use default [[
    subset_hhh4sims_attributes(x, i, j)
}

## select a specific "hhh4sims" from the list of simulations
## (the inverse of as.hhh4simslist.hhh4sims(xx))
"[[.hhh4simslist" <- function (x, i)
{
    xx <- NextMethod("[[")
    a <- attributes(xx)
    attributes(xx) <- c(a[c("dim", "dimnames")],
                        attributes(x)[c("initial", "stsObserved")],
                        list(class = "hhh4sims"),
                        a[c("call", "seed")])
    xx
}

## aggregate predictions over time and/or (groups of) units
aggregate.hhh4simslist <- function (x, units = TRUE, time = FALSE, ..., drop = FALSE)
{
    if (drop || time) { # unclass(x) to use default "[["-method in lapply
        lapply(X = unclass(x), FUN = aggregate.hhh4sims,
               units = units, time = time, ..., drop = TRUE)
    } else {
        as.hhh4simslist.list(
            lapply(X = x, FUN = aggregate.hhh4sims,
                   units = units, time = time, ..., drop = FALSE)
            )
    }
}


####################
### plot methods ###
####################

check_groups <- function (groups, units)
{
    if (is.null(groups)) {
        factor(rep.int("overall", length(units)))
    } else if (isTRUE(groups)) {
        factor(units, levels = units)
    } else {
        stopifnot(length(groups) == length(units))
        as.factor(groups)
    }
}

plot.hhh4simslist <- function (x, type = c("size", "time", "fan"), ...,
                               groups = NULL, par.settings = list())
{
    FUN <- paste("plotHHH4sims", match.arg(type), sep = "_")
    groups <- check_groups(groups, colnames(attr(x, "stsObserved"), do.NULL=FALSE))
    ngroups <- nlevels(groups)
    if (is.list(par.settings)) {
        par.defaults <- list(mar = c(4,4,2,0.5)+.1, las = 1)
        if (ngroups > 1)
            par.defaults$mfrow <- sort(n2mfrow(ngroups))
        par.settings <- modifyList(par.defaults, par.settings)
        opar <- do.call("par", par.settings)
        on.exit(par(opar))
    }

    if (ngroups == 1) {
        do.call(FUN, list(quote(x), ...))
    } else { # stratified plots by groups of units
        invisible(sapply(
            X = levels(groups),
            FUN = function (group) {
                x_group <- x[, which(group == groups) , ] # [-method has drop=F
                do.call(FUN, list(quote(x_group), ..., main = group))
            },
            simplify = FALSE, USE.NAMES = TRUE))
    }
}


### simulated final size distribution as boxplots aggregated over all units

plotHHH4sims_size <- function (x, horizontal = TRUE, trafo = NULL,
                               observed = TRUE, names = base::names(x), ...)
{
    x <- as.hhh4simslist(x)
    if (horizontal) {
        names <- rev(names)
        x <- rev(x)
    }
    if (is.null(trafo)) #trafo <- scales::identity_trans()
        trafo <- list(name = "identity", transform = identity)
    if (isTRUE(observed)) observed <- list()
    nsims <- sapply(X = unclass(x), # simply use the default "[["-method
                    FUN = colSums, dims = 2, # sum over 1:2 (time x unit)
                    simplify = TRUE, USE.NAMES = TRUE)
    nsimstrafo <- trafo$transform(nsims)

    ## default boxplot arguments
    fslab <- "size"
    if (trafo$name != "identity")
        fslab <- paste0(fslab, " (", trafo$name, "-scale)")
    defaultArgs <- list(ylab=fslab, yaxt="n", las=1, cex.axis=1, border=1)
    if (horizontal) names(defaultArgs) <- sub("^y", "x", names(defaultArgs))
    ## defaultArgs$mai <- par("mai")
    ## defaultArgs$mai[2] <- max(strwidth(boxplot.args$names, units="inches",
    ##                                    cex=boxplot.args$cex.axis))
    ## if (trafo$name != "identity") {
    ##     ## ?bxp: 'yaxs' and 'ylim' are used 'along the boxplot'
    ##     defaultArgs <- c(defaultArgs,
    ##                      list(ylim=c(0,max(nsimstrafo)*1.05), yaxs="i"))
    ## }

    ## generate boxplots
    boxplot.args <- modifyList(defaultArgs, list(...))
    boxplot.args$horizontal <- horizontal
    boxplot.args$names <- names
    do.call("boxplot", c(list(x=nsimstrafo), boxplot.args))

    ## add means
    if (horizontal) {
        points(x=colMeans(nsimstrafo), y=1:ncol(nsimstrafo), pch=8, col=boxplot.args$border)
    } else points(colMeans(nsimstrafo), pch=8, col=boxplot.args$border)

    ## add axis
    aty <- pretty(nsims, n=par("lab")[2-horizontal])
    ##aty <- checkat(list(n=par("lab")[2], trafo=trafo), nsims) # linear on sqrt-scale
    axis(2-horizontal, at=trafo$transform(aty), labels=aty, las=boxplot.args$las)

    ## add line showing observed size
    if (is.list(observed)) {
        nObs <- sum(observed(attr(x, "stsObserved")))
        observed <- modifyList(
            list(col = 1, lty = 2, lwd = 2,
                 labels = nObs, font = 2, las = boxplot.args$las,
                 mgp = if (horizontal) c(3, 0.4, 0)),
            observed)
        observed_line <- c(
            setNames(list(trafo$transform(nObs)), if (horizontal) "v" else "h"),
            observed[c("col", "lty", "lwd")])
        do.call("abline", observed_line)
        if (!is.null(observed[["labels"]]))
            do.call("axis", c(
                list(side = 2-horizontal, at = trafo$transform(nObs)),
                observed))
    }

    ## numeric summary
    mysummary <- function(x)
        c(mean=mean(x), quantile(x, probs=c(0.025, 0.5, 0.975)))
    nsum <- t(apply(nsims, 2, mysummary))
    invisible(nsum)
}


### Plot mean time series of the simulated counts

plotHHH4sims_time <- function (
    x, average = mean, individual = length(x) == 1,
    conf.level = if (individual) 0.95 else NULL, #score = "rps",
    matplot.args = list(), initial.args = list(), legend = length(x) > 1,
    xlim = NULL, ylim = NULL, add = FALSE, ...)
{
    x <- as.hhh4simslist(x)
    nModels <- length(x)
    ytInit <- rowSums(attr(x, "initial"))
    stsObserved <- attr(x, "stsObserved")
    ytObs <- rowSums(observed(stsObserved))
    ytSim <- aggregate.hhh4simslist(x, units = TRUE, time = FALSE, drop = TRUE)
    average <- match.fun(average)
    ytMeans <- vapply(
        X = ytSim,
        FUN = function (x) apply(x, 1, average),
        FUN.VALUE = numeric(length(ytObs)), USE.NAMES = TRUE)

    ## axis range
    if (is.null(xlim) && is.list(initial.args))
        xlim <- c(1 - length(ytInit) - 0.5, length(ytObs) + 0.5)
    if (is.null(ylim))
        ylim <- c(0, max(ytObs, if (individual)
            unlist(ytSim, recursive = FALSE, use.names = FALSE) else ytMeans))

    ## graphical parameters
    stopifnot(is.list(matplot.args))
    matplot.args <- modifyList(
        list(y = ytMeans, type = "b", lty = 1, lwd = 3, pch = 20,
             col = rainbow(nModels)),
        matplot.args)
    col <- rep_len(matplot.args$col, nModels)

    ## observed time series data during simulation period
    if (!add)
        plot(stsObserved, type = observed ~ time, xlim = xlim, ylim = ylim, ...)

    ## add initial counts
    if (is.list(initial.args)) {
        initial.args <- modifyList(
            list(x = seq(to = 0, by = 1, length.out = length(ytInit)),
                 y = ytInit, type = "h", lwd = 5),
            initial.args)
        do.call("lines", initial.args)
    }

    ## add counts of individual simulation runs
    if (individual) {
        for (i in seq_len(nModels))
            matlines(ytSim[[i]], lty=1, col=if (requireNamespace("scales"))
                scales::alpha(col[i], alpha=0.1) else col[i])
        col <- col2rgb(col)
        col <- apply(col, 2, function (x)
                     if (all(x == 0)) "grey" else
                     do.call("rgb", as.list(x / 255 * 0.5)))
    }

    ## add means (or medians)
    matplot.args[["col"]] <- col
    do.call("matlines", matplot.args)

    ## add CIs
    if (isScalar(conf.level)) {
        alpha2 <- (1-conf.level)/2
        ytQuant <- lapply(ytSim, function (sims)
            t(apply(sims, 1, quantile, probs=c(alpha2, 1-alpha2))))
        matlines(sapply(ytQuant, "[", TRUE, 1L), col=col, lwd=matplot.args$lwd, lty=2)
        matlines(sapply(ytQuant, "[", TRUE, 2L), col=col, lwd=matplot.args$lwd, lty=2)
    }

    ## add scores
    ## if (length(score)==1) {
    ##     scorestime <- simplify2array(
    ##         simscores(x, by="time", scores=score, plot=FALSE),
    ##         higher=FALSE)
    ##     matlines(scales::rescale(scorestime, to=ylim),
    ##              lty=2, lwd=1, col=col)
    ## }

    ## add legend
    if (!identical(FALSE, legend)) {
        xnames <- if (is.vector(legend, mode = "character")) {
            if (length(legend) != length(x))
                warning("'length(legend)' should be ", length(x))
            legend
        } else {
            names(x)
        }
        legendArgs <- list(x="topright", legend=xnames, bty="n",
                           col=col, lwd=matplot.args$lwd, lty=matplot.args$lty)
        if (is.list(legend))
            legendArgs <- modifyList(legendArgs, legend)
        do.call("legend", legendArgs)
    }

    ## Done
    ret <- cbind(observed = ytObs, ytMeans)
    ## if (length(score) == 1)
    ##     attr(ret, score) <- scorestime
    invisible(ret)
}


### Better for a single model: "fanplot"

plotHHH4sims_fan <- function (x, which = 1,
    fan.args = list(), observed.args = list(), initial.args = list(),
    means.args = NULL, key.args = NULL, xlim = NULL, ylim = NULL,
    add = FALSE, xaxis = list(), ...)
{
    x <- as.hhh4simslist(x)[[which]]
    ytInit <- rowSums(attr(x, "initial"))
    stsObserved <- attr(x, "stsObserved")
    ytObs <- rowSums(observed(stsObserved))
    ytSim <- aggregate.hhh4sims(x, units = TRUE, time = FALSE, drop = TRUE)

    ## graphical parameters
    if (is.null(xlim) && is.list(initial.args))
        xlim <- c(1 - length(ytInit) - 0.5, length(ytObs) + 0.5)
    stopifnot(is.list(fan.args))
    fan.args <- modifyList(
        list(probs = seq.int(0.01, 0.99, 0.01)),
        fan.args, keep.null = TRUE)

    ## compute the quantiles
    quantiles <- t(apply(ytSim, 1, quantile, probs = fan.args$probs))

    ## create (or add) the fanplot
    fanplot(quantiles = quantiles, probs = fan.args$probs,
            means = rowMeans(ytSim), observed = ytObs,
            fan.args = fan.args, means.args = means.args,
            observed.args = observed.args, key.args = key.args,
            xlim = xlim, ylim = ylim, add = add,
            xaxt = if (is.list(xaxis)) "n" else "s", ...)

    ## add initial counts
    if (is.list(initial.args)) {
        initial.args <- modifyList(
            list(x = seq(to = 0, by = 1, length.out = length(ytInit)),
                 y = ytInit, type = "p", pch = 19),
            initial.args)
        do.call("lines", initial.args)
    }

    ## add time axis
    if (is.list(xaxis)) {
        xaxis <- modifyList(list(epochsAsDate = TRUE), xaxis)
        do.call("addFormattedXAxis", c(list(x = stsObserved), xaxis))
    }

    invisible(NULL)
}
