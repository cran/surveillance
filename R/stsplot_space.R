################################################################################
### Snapshot map (spplot) of an sts-object or matrix of counts
###
### Copyright (C) 2013-2014,2016,2017,2020,2021 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################

## x: "sts" or (simulated) matrix of counts
## tps: one or more time points. The unit-specific _sum_ of time points "tps" is
##      plotted. tps=NULL means cumulation over all time points in x.
## at: number of levels for the grouped counts or specific break points to
##     use, or list(n, data, trafo) passed to getPrettyIntervals(),
##     where data and trafo are optional.
##     CAVE: intervals are closed on the left and open to the right.
##     From panel.levelplot: zcol[z >= at[i] & z < at[i + 1]] <- i
##     i.e. at=0:1 will have NA (=also white) for counts=1, thus we have to
##     ensure max(at) > max(counts)

stsplot_space <- function (x, tps = NULL, map = x@map, population = NULL,
                           main = NULL, labels = FALSE,
                           ..., # placed here to support passing 'col'
                           at = 10, col.regions = NULL,
                           colorkey = list(space="bottom", labels=list(at=at)),
                           total.args = NULL,
                           gpar.missing = list(col="darkgrey", lty=2, lwd=2),
                           sp.layout = NULL,
                           ## aspect = mapasp(map), # currently hard-coded
                           xlim = bbox(map)[1, ], ylim = bbox(map)[2, ])
{
    counts <- if (inherits(x, "sts")) observed(x) else x
    if (is.null(tps))
        tps <- seq_len(nrow(counts))
    if (length(map) == 0L) stop("no map")
    if (is.null(colnames(counts)))
        stop("need 'colnames(x)' (to be matched against 'row.names(map)')")
    if (!all(colnames(counts) %in% row.names(map)))
        stop("incomplete 'map'; ensure that 'all(colnames(x) %in% row.names(map))'")

    ## compute data to plot
    ncases <- getCumCounts(counts, tps)
    total <- sum(ncases)
    if (!is.null(population)) { # divide counts by region-specific population
        population <- parse_population_argument(population, x)  # pop matrix
        populationByRegion <- population[tps[1L],] # pop at first time point
        ncases <- ncases / populationByRegion  # (cumulative) incidence by region
        total <- total / sum(populationByRegion)
    }

    ## add ncases to map@data
    map <- as(map, "SpatialPolygonsDataFrame")
    map$ncases <- NA_real_
    map$ncases[match(colnames(counts),row.names(map))] <- ncases

    ## default main title
    if (is.null(main) && inherits(x, "sts"))
        main <- stsTimeRange2text(x, tps)

    ## check/determine color break points 'at'
    at <- checkat(at, ncases, counts = is.null(population))
    ## default color palette
    if (is.null(col.regions)) {
        separate0 <- is.null(population) && at[1] == 0 && at[2] <= 1
        col.regions <- c(if (separate0) "white",
                         .hcl.colors(length(at)-1-separate0))
    }
    ## colorkey settings
    if (!missing(colorkey) && is.list(colorkey))
        colorkey <- modifyList(eval(formals()$colorkey), colorkey)

    ## automatic additions to sp.layout (region labels and total)
    if (is.list(gpar.missing) && anyNA(map$ncases)) {
        layout.missing <- c(list("sp.polygons", obj=map[is.na(map$ncases),]),
                            gpar.missing)
        sp.layout <- c(sp.layout, list(layout.missing))
    }
    if (!is.null(layout.labels <- layout.labels(map, labels))) {
        sp.layout <- c(sp.layout, list(layout.labels))
    }
    if (is.list(total.args)) {
        total.args <- modifyList(list(label="Overall: ", x=1, y=0),
                                 total.args)
        if (is.null(total.args$just))
            total.args$just <- with (total.args, if (all(c(x,y) %in% 0:1)) {
                c(c("left", "right")[1+x], c("bottom","top")[1+y])
            } else "center")
        total.args$label <- paste0(total.args$label, round(total,1))
        layout.total <- c(grid::grid.text, total.args)
        ## "grid.text" wouldn't work since package "sp" doesn't import it
        sp.layout <- c(sp.layout, list(layout.total))
    }

    ## generate the spplot()
    args <- list(quote(map[!is.na(map$ncases),]), "ncases", main=main,
                 col.regions=col.regions, at=at, colorkey=colorkey,
                 sp.layout=sp.layout, aspect=mapasp(map), xlim=xlim, ylim=ylim,
                 quote(...))
    do.call("spplot", args)
}



#######################################################
### Auxiliary functions for the "sts" snapshot function
#######################################################

## sum of counts by unit over time points "tps"
## the resulting vector has no names
getCumCounts <- function (counts, tps)
{
    ntps <- length(tps)
    if (ntps == 1) {
        counts[tps,]
    } else {
        .colSums(counts[tps,,drop=FALSE], ntps, ncol(counts))
    }
}

parse_population_argument <- function (population, x)
{
    if (is.matrix(population)) {
        if (!identical(dim(population), dim(x)))
            stop("'dim(population)' does not match the data dimensions")
    } else if (isScalar(population)) { # a unit, e.g., per 1000 inhabitants
        if (!inherits(x, "sts"))
            stop("'", deparse(substitute(x)), "' is no \"sts\" object; ",
                 "population numbers must be supplied")
        population <- population(x) / population
    } else { # region-specific population numbers (as in surveillance <= 1.12.2)
        stopifnot(is.vector(population, mode = "numeric"))
        if (length(population) != ncol(x))
            stop("'length(population)' does not match the number of data columns")
        population <- rep(population, each = nrow(x))
        dim(population) <- dim(x)
    }
    population
}

checkat <- function (at, data, counts = TRUE) { # for non-transformed "data"
    if (isScalar(at))
        at <- list(n=at)
    if (is.list(at)) {
        at <- modifyList(list(n=10, data=data, counts=counts), at)
        do.call("getPrettyIntervals", at)
    } else { # manual breaks
        stopifnot(is.vector(at, mode = "numeric"), !anyNA(at))
        at <- sort(at)
        r <- range(data, na.rm = TRUE)
        c(if (r[1L] < at[1L]) 0,
          at,
          if (r[2L] >= at[length(at)]) {
              ## round up max to 1 significant digit (including 0.1 to 0.2)
              .decs <- 10^floor(log10(r[2L]))
              ceiling(r[2L]/.decs + sqrt(.Machine$double.eps))*.decs
          })
    }
}

getPrettyIntervals <- function (n, data, trafo=NULL, counts=TRUE, ...) {
    maxcount <- max(data, na.rm=TRUE)
    if (counts && maxcount < n) { # no aggregation of counts necessary
        at <- 0:ceiling(maxcount+sqrt(.Machine$double.eps)) # max(at) > maxcount
    } else {
        at <- if (is.null(trafo)) { # equivalent to trafo=scales::sqrt_trans()
            pretty(sqrt(data), n=n+1, ...)^2
        } else {
            scales::trans_breaks(trafo$trans, trafo$inv, n=n+1, ...)(data)
        }
        ## { # alternative: quantile-based scale (esp. for incidence plots)
        ##     quantile(data, probs=seq(0,1,length.out=n+1), na.rm=TRUE)
        ## }
        if (counts && at[1] == 0 && at[2] > 1) # we want 0 counts separately ("white")
            at <- sort(c(1, at))
        if (at[length(at)] == maxcount) # ensure max(at) > max(data)
            at[length(at)] <- at[length(at)] + if (counts) 1 else 0.001*diff(range(at))
    }
    at
}

stsTime2text <- function (stsObj, tps=TRUE, fmt=NULL)
{
    if (is.null(fmt))
        fmt <- switch(as.character(stsObj@freq),
                      "1" = "%i", "52" = "%i-W%02i", "%i/%i")
    sprintf(fmt, year(stsObj)[tps], epochInYear(stsObj)[tps])
}

stsTimeRange2text <- function (stsObj, tps, fmt=NULL, sep=" to ")
{
    tpsRangeYW <- stsTime2text(stsObj, tps=range(tps), fmt=fmt)
    paste0(unique(tpsRangeYW), collapse=sep)
}
