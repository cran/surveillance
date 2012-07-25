################################################################################
### Function 'simEpidataCS' simulates a point pattern according to an
### additive-multiplicative spatio-temporal intensity model of class 'twinstim'.
### It basically uses Ogata's modified thinning algorithm
### (cf. Daley & Vere-Jones, 2003, Algorithm 7.5.V.).
### Author: Sebastian Meyer
################################################################################

# CAVE: the type of contrasts for factor variables has to be set through options("contrasts")

### TODO: if epidemic-only process (!hash), we actually don't need stgrid
###       the function is not yet arranged for endemic-only or epidemic-only settings

simEpidataCS <- function (endemic, epidemic, siaf, tiaf, qmatrix, rmarks,
    events, stgrid, tiles, beta0, beta, gamma, siafpars, tiafpars,
    t0 = stgrid$start[1], T = tail(stgrid$stop,1), nEvents = 1e5, nCub,
    W = NULL, trace = 5, nCircle2Poly = 32, gmax = NULL, .allocate = 500,
    .skipChecks = FALSE, .onlyEvents = FALSE, ...)
{
    ptm <- proc.time()[[3]]
    cl <- match.call()


    #######################
    ### Check arguments ### (this takes many lines of code ...)
    #######################


    cat("\nChecking the supplied arguments ...\n")

    ### Some simple checks

    if (length(trace) != 1L) stop("'trace' must be a single integer or logical value")
    trace <- as.integer(trace)
    if (!isScalar(nCircle2Poly)) stop("'nCircle2Poly' must be scalar")
    nCircle2Poly <- as.integer(nCircle2Poly)
    if (!isScalar(.allocate)) stop("'.allocate' must be scalar")
    .allocate <- as.integer(.allocate)
    .skipChecks <- as.logical(.skipChecks)
    .onlyEvents <- as.logical(.onlyEvents)


    ### Check formulae

    if (missing(endemic)) endemic <- ~ 0 else stopifnot(inherits(endemic, "formula"))
    if (missing(epidemic)) epidemic <- ~ 0 else stopifnot(inherits(epidemic, "formula"))


    ### Check parameters

    beta0 <- if (missing(beta0)) numeric(0L) else as.vector(beta0, mode="numeric")
    beta  <- if (missing(beta))  numeric(0L) else as.vector(beta,  mode="numeric")
    gamma <- if (missing(gamma)) numeric(0L) else as.vector(gamma, mode="numeric")
    siafpars <- if (missing(siafpars)) numeric(0L) else as.vector(siafpars, mode="numeric")
    tiafpars <- if (missing(tiafpars)) numeric(0L) else as.vector(tiafpars, mode="numeric")
    nbeta0 <- length(beta0)
    p <- length(beta)
    q <- length(gamma)
    nsiafpars <- length(siafpars)
    ntiafpars <- length(tiafpars)

    hase <- q > 0L
    hassiafpars <- nsiafpars > 0L
    hastiafpars <- ntiafpars > 0L
    if (!hase && (hassiafpars | hastiafpars)) {
        stop("'siafpars' and 'tiafpars' require 'gamma'")
    }


    ### Check qmatrix

    if (missing(qmatrix)) qmatrix <- diag(1)
    typeNames <- rownames(qmatrix)
    if (is.null(typeNames)) typeNames <- "1"    # single event type
    qmatrix <- checkQ(qmatrix, typeNames)
    nTypes <- length(typeNames)
    if (nbeta0 > 1L && nbeta0 != nTypes) {
        stop("'beta0' must have length 0, 1, or 'nrow(qmatrix)'")
    }
    qSumTypes <- rowSums(qmatrix)  # how many types can be triggered by each type


    ### Check stgrid

    if (!.skipChecks) {
        cat("Checking 'stgrid':\n")
        stgrid <- checkstgrid(stgrid[grep("^BLOCK$", names(stgrid), invert=TRUE)])
    }


    ### Check time range

    if (is.null(t0)) t0 <- eval(formals()$t0)
    if (is.null(T)) T <- eval(formals()$T)
    if (!isScalar(t0) || !isScalar(T)) {
        stop("endpoints 't0' and 'T' must be single numbers")
    }
    if (T <= t0) {
        stop("'T' must be greater than 't0'")
    }
    stopifnot(t0 >= stgrid$start[1], T <= tail(stgrid$stop,1))


    ### Subset stgrid to include actual time range only

    # BLOCK in stgrid such that start time is equal to or just before t0
    block_t0 <- stgrid$BLOCK[match(TRUE, stgrid$start > t0) - 1L]
    # BLOCK in stgrid such that stop time is equal to or just after T
    block_T <- stgrid$BLOCK[match(TRUE, stgrid$stop >= T)]
    stgrid <- subset(stgrid, BLOCK >= block_t0 & BLOCK <= block_T)
    stgrid$start[stgrid$BLOCK == block_t0] <- t0
    stgrid$stop[stgrid$BLOCK == block_T] <- T
    # matrix of BLOCKS and start times (used later)
    blockstarts <- with(stgrid,
        cbind(block_t0:block_T, start[match(block_t0:block_T, BLOCK)])
    )


    ### Check class of W, and class, proj4string and names of tiles

    if (is.null(W)) {
    	if (require("maptools")) {
            cat("Building W as the union of 'tiles' ...\n")
            W <- maptools::unionSpatialPolygons(tiles,
                     IDs = rep.int(1,length(tiles@polygons)),
                     avoidGEOS = TRUE)
        } else {
            stop("automatic generation of 'W' from 'tiles' requires package \"maptools\"")
        }
    }
    stopifnot(inherits(W, "SpatialPolygons"), inherits(tiles, "SpatialPolygons"),
              proj4string(tiles) == proj4string(W),
              (tileLevels <- levels(stgrid$tile)) %in% row.names(tiles))

    # Transform W into a gpc.poly
    Wgpc <- as(W, "gpc.poly")


    ### Check mark-generating function

    epscols <- c("eps.t", "eps.s")
    unpredMarks <- setdiff(unique(c(if (hase) epscols, all.vars(epidemic))), "type")
    # <- eps.t and eps.s are also taken to be unpredictable marks
    rmarks <- if (missing(rmarks) || is.null(rmarks)) NULL else match.fun(rmarks)
    hasMarks <- !is.null(rmarks)
    if (hasMarks) {
        sampleCoordinate <- coordinates(spsample(tiles, n=1L, type="random"))
        sampleMarks <- rmarks(t0, sampleCoordinate)
        # should be a one-row data.frame
        if (!is.data.frame(sampleMarks) || nrow(sampleMarks) != 1L) {
            stop("'rmarks' must return a one-row data.frame of marks")
        }
        if (!all(sapply(sampleMarks, function(x) inherits(x, c("integer","numeric","factor"), which=FALSE)))) {
            stop("'rmarks' must return \"numeric\", \"integer\" or \"factor\" variables only")
        }
        markNames <- names(sampleMarks)
        if (.idx <- match(FALSE, unpredMarks %in% markNames, nomatch=0L)) {
            stop("the unpredictable mark '", unpredMarks[.idx], "' is not returned by 'rmarks'")
        }
    } else if (length(unpredMarks) > 0L) {
        stop("missing mark-generating function 'rmarks' for the 'epidemic' component")
    }


    ### Check events (prehistory of the process)

    # empty prehistory
    eventCoords <- matrix(0, nrow=0L, ncol=2L)
    eventData <- data.frame(
        ID   = integer(0L),
        time = numeric(0L),
        tile = factor(character(0L), levels=tileLevels),
        type = factor(character(0L), levels=typeNames),
        check.rows = FALSE, check.names = FALSE
    )
    if (hasMarks) {
        eventData <- cbind(eventData, sampleMarks[0L,])
    }
    Nout <- 0L
    # do we have a prehistory in 'events'
    if (!missing(events) && !is.null(events)) {
        .reservedColsIdx <- match(c("ID",".obsInfLength",".bdist",".influenceRegion",".sources"), names(events))
        events <- events[setdiff(seq_along(names(events)),.reservedColsIdx)]
        if (!.skipChecks) {
            cat("Checking 'events':\n")
            events <- checkEvents(events, dropTypes = FALSE)
            # epscols are obligatory in 'checkEvents', which is also appropriate here
            # because in the case of endemic-only processes we would not need a prehistory at all
        }
        # select prehistory of events which are still infective
        .stillInfective <- with(events@data, time <= t0 & time + eps.t > t0)
        Nout <- sum(.stillInfective)    # = number of events in the prehistory
        if (Nout > 0L) {
            stopifnot(proj4string(events) == proj4string(W))
            events <- events[.stillInfective,]
            # check event types
            events@data$type <- factor(events@data$type, levels=typeNames)
            if (any(.typeIsNA <- is.na(events@data$type))) {
                warning("removed unknown event types in 'events'")
                events <- events[!.typeIsNA,]
            }
            eventCoords <- coordinates(events)
            eventData <- events@data
            # update ID counter
            eventData$ID <- 1:Nout
            # check presence of unpredictable marks
            if (.idx <- match(FALSE, unpredMarks %in% names(eventData), nomatch=0L)) {
                stop("missing unpredictable mark '", unpredMarks[.idx], "' in 'events'")
            }
            # check type of unpredictable marks
            for (um in unpredMarks) {
                if (!identical(class(sampleMarks[[um]]), class(eventData[[um]]))) {
                    stop("the class of the unpredictable mark '", um, "' in the 'events' prehistory ",
                         "is not identical to the class returned by 'rmarks'")
                }
            }
            # add marks which are not in epidemic model but which are simulated by 'rmarks'
            if (length(.add2events <- setdiff(markNames, names(eventData))) > 0L) {
                eventData <- cbind(eventData, sampleMarks[.add2events])
                is.na(eventData[.add2events]) <- TRUE
            }
            eventData <- eventData[c("ID", "time", "tile", "type", markNames)]
        } else {
            .eventstxt <- if (.skipChecks) "data$events" else "events"   # account for simulate.twinstim
            cat("(no events from '", .eventstxt, "' were considered as prehistory)\n", sep="")
        }
    }


    ### Build epidemic model matrix

    epidemic <- terms(epidemic, data = eventData, keep.order = TRUE)
    if (!is.null(attr(epidemic, "offset"))) {
        warning("offsets are not implemented for the 'epidemic' component")
    }

    # helper function taking eventData and returning the epidemic model.matrix
    buildmme <- function (eventData)
    {
        mfe <- model.frame(epidemic, data = eventData, na.action = na.fail, drop.unused.levels = FALSE)
        model.matrix(epidemic, mfe)
    }
    mme <- buildmme(eventData)

    if (ncol(mme) != q) {
        cat(ncol(mme), "epidemic model terms:\t", paste(colnames(mme), collapse="  "), "\n")
        stop("length of 'gamma' (", q, ") does not match the 'epidemic' specification (", ncol(mme), ")")
    }


    ### Build endemic model matrix

    endemic <- terms(endemic, data = stgrid, keep.order = TRUE)

    # check if we have an endemic component at all
    hasOffset <- !is.null(attr(endemic, "offset"))
    hash <- (nbeta0 + p + hasOffset) > 0L
    if (!hash && !hase) {
        stop("nothing to do: neither endemic nor epidemic parameters were specified")
        # actually, the process might be endemic offset-only, which I don't care about ATM
    }

    # remove (1|type) specification
    typeSpecificEndemicIntercept <-
        "1 | type" %in% attr(endemic, "term.labels") || nbeta0 > 1
    if (typeSpecificEndemicIntercept) {
        endemic <- update(endemic, ~ . - (1|type)) # this drops the terms attributes
        endemic <- terms(endemic, data = data$stgrid, keep.order = TRUE)
        if (nbeta0 < 1L) {
            stop("for type-specific endemic intercepts, 'beta0' must be longer than 1")
        }
    }

    # ensure that we have correct contrasts in the endemic component
    attr(endemic, "intercept") <- as.integer(nbeta0 > 0L)

    # which variables do we have to copy from stgrid?
    stgridCopyCols <- match(all.vars(endemic), names(stgrid), nomatch = 0L)

    # helper function taking eventData (with time and tile columns) and returning the endemic model.matrix
    buildmmh <- function (eventData)
    {
        # if 'pi' appears in 'endemic' we don't care, and if a true covariate is missing, model.frame will throw an error
        # attach endemic covariates from 'stgrid' to events
        gridcellsOfEvents <- integer(nrow(eventData))
        for (i in seq_along(gridcellsOfEvents)) {
            gridcellsOfEvents[i] <- gridcellOfEvent(eventData[i,"time"], eventData[i,"tile"], stgrid)
        }
        eventData <- cbind(eventData, stgrid[gridcellsOfEvents, stgridCopyCols])
        # construct model matrix
        mfhEvents <- model.frame(endemic, data = eventData, na.action = na.fail, drop.unused.levels = FALSE)
        mmhEvents <- model.matrix(endemic, mfhEvents)
        # exclude intercept from endemic model matrix below, will be treated separately
        if (nbeta0 > 0) mmhEvents <- mmhEvents[,-1,drop=FALSE]
        structure(mmhEvents, offset = model.offset(mfhEvents))
    }

    # actually, we don't need the endemic model matrix for the events at all
    # this is just to test consistence with 'beta'
    mmh <- buildmmh(eventData[0L,])
    if (ncol(mmh) != p) {
        stop("length of 'beta' (", p, ") does not match the 'endemic' specification (", ncol(mmh), ")")
    }


    ### Build endemic model matrix on stgrid

    mfhGrid <- model.frame(endemic, data = stgrid, na.action = na.fail, drop.unused.levels = FALSE,
                           BLOCK = BLOCK, tile = tile, ds = area)
    # we don't actually need 'tile' in mfhGrid; this is only for easier identification when debugging
    mmhGrid <- model.matrix(endemic, mfhGrid)
    # exclude intercept from endemic model matrix below, will be treated separately
    if (nbeta0 > 0) mmhGrid <- mmhGrid[,-1,drop=FALSE]

    # Extract endemic model components
    offsetGrid <- model.offset(mfhGrid)
    gridBlocks <- mfhGrid[["(BLOCK)"]]
    ds <- mfhGrid[["(ds)"]]


    ### Parse interaction functions

    if (hase) {

        ## Check interaction functions
        siaf <- do.call(".parseiaf", args = alist(siaf))
        tiaf <- do.call(".parseiaf", args = alist(tiaf))

        ## Spatially constant interaction siaf(s) = 1
        constantsiaf <- is.null(siaf) || attr(siaf, "constant")
        if (constantsiaf) siaf <- siaf.constant()
        if (siaf$npars != nsiafpars) {
            stop("length of 'siafpars' (", nsiafpars,
                 ") does not match the 'siaf' specification (", siaf$npars, ")")
        }

        ## Temporally constant interaction tiaf(t) = 1
        constanttiaf <- is.null(tiaf) || attr(tiaf, "constant")
        if (constanttiaf) {
            tiaf <- tiaf.constant()
            gmax <- 1L
        }
        if (tiaf$npars != ntiafpars) {
            stop("length of 'tiafpars' (", ntiafpars,
                 ") does not match the 'tiaf' specification (", tiaf$npars, ")")
        }

        ## Check nCub
        if (!constantsiaf) {
            nCub <- as.integer(nCub)
            if (any(nCub <= 0L)) {
                stop("'nCub' must be positive")
            }
        }

        ## Check gmax
        if (is.null(gmax)) {
            gmax <- max(tiaf$g(rep.int(0,nTypes), tiafpars, 1:nTypes))
            cat("assuming gmax =", gmax, "\n")
        } else if (!isScalar(gmax)) {
            stop("'gmax' must be scalar")
        }

    } else {
        if ((!missing(siaf) && !is.null(siaf)) || (!missing(tiaf) && !is.null(tiaf))) {
            warning("'siaf' and 'tiaf' can only be modelled in conjunction with an 'epidemic' process")
        }
        siaf <- tiaf <- NULL
    }


    ### print some information on the upcoming simulation

    txtPrehistory <- if (Nout == 0L) "no prehistory" else paste(Nout,
        ngettext(Nout, "event", "events"), "in the prehistory")
    cat("\nSimulating a", if (hasMarks) "marked", "spatio-temporal point pattern with",
        "\n\t-", nTypes, ngettext(nTypes, "event type", "event types"),
        "\n\t-", txtPrehistory)

    coefs <- c(
        if (nbeta0 > 1L) {
            structure(beta0, names=paste0("h.type",typeNames))
        } else if (nbeta0 == 1L) structure(beta0, names="h.(Intercept)"),
        if (p > 0L) structure(beta, names=paste("h",colnames(mmh),sep=".")),
        if (hase) structure(gamma, names=paste("e",colnames(mme),sep=".")),
        if (hassiafpars) structure(siafpars, names=paste("e.siaf",1:nsiafpars,sep=".")),
        if (hastiafpars) structure(tiafpars, names=paste("e.tiaf",1:ntiafpars,sep="."))
    )
    cat("\n\t-", length(coefs), "coefficients:\n\n")
    print(coefs)



    ##########################################
    ### CIF of the temporal ground process ###
    ##########################################


    ### calculate integral of endemic component over W (= union of tiles) and types for all time blocks in stgrid

    dsexpeta <- local({
        eta <- drop(mmhGrid %*% beta)
        if (!is.null(offsetGrid)) eta <- offsetGrid + eta
        ds * exp(unname(eta))
    })
    fact <- if (nbeta0 > 1L) sum(exp(beta0)) else if (nbeta0 == 1L) nTypes*exp(unname(beta0)) else nTypes
    hIntWK <- fact * c(tapply(dsexpeta, gridBlocks, sum))
    #<- is a named vector with names referencing BLOCK in stgrid


    ### calculate terms appearing in the epidemic component (the lambdag function uses eTerms)

    # compute computationally effective range of the 'siaf' function
    if (!constantsiaf && !is.null(siaf$Fcircle)) {
        effRangeTypes <- rep(siaf$effRange(siafpars), length.out = nTypes)
    }

    # helper function evaluating the epidemic terms of the ground intensity for a specific set of events
    eTermsCalc <- function (eventData, eventCoords)
    {
        # epidemic model matrix (will be multiplied with gamma)
        mme <- buildmme(eventData)
        # distance to the border (required for siafInt below, and for epidataCS)
        .bdist <- bdist(eventCoords, Wgpc)
        # spatial influence regions of the events
        .influenceRegion <- if (nrow(eventCoords) > 0L) .influenceRegions(
            events = SpatialPointsDataFrame(
                coords = eventCoords,
                data = data.frame(eps.s = eventData$eps.s, .bdist = .bdist),
                match.ID = FALSE
            ),
            Wgpc = Wgpc,
            npoly = nCircle2Poly
        ) else list()
        # integrate the two-dimensional 'siaf' function over the influence region
        siafInt <- if (length(.influenceRegion) == 0L) numeric(0L) else if (constantsiaf) {
                sapply(.influenceRegion, attr, "area")
            } else if (is.null(siaf$Fcircle)) { # if siaf$Fcircle is not available
                sapply(seq_along(.influenceRegion), function (i) {
                    polyCub.midpoint(.influenceRegion[[i]], siaf$f, siafpars, as.integer(eventData[i,"type"]), eps = nCub)
                })
            } else { # fast integration over circular domains !
                # effRange of each event
                effRanges <- effRangeTypes[eventData$type]
                # automatic choice of h (pixel spacing in image)
                hs <- effRanges / nCub   # could e.g. equal sigma
                # Compute the integral of 'siaf' over each influence region
                siafInts <- numeric(length(.influenceRegion))
                for(i in seq_along(siafInts)) {
                    eps <- eventData$eps.s[i]
                    bdisti <- .bdist[i]
                    effRange <- effRanges[i]
                    siafInts[i] <- if (effRange <= bdisti) {
                            # effective region ("6 sigma") completely inside W
                            siaf$Fcircle(min(eps,effRange), siafpars, as.integer(eventData[i,"type"]))
                        } else if (eps <= bdisti) {
                            # influence region is completely inside W
                            siaf$Fcircle(eps, siafpars, as.integer(eventData[i,"type"]))
                        } else {
                            # integrate over polygonal influence region
                            polyCub.midpoint(.influenceRegion[[i]], siaf$f,
                                siafpars, as.integer(eventData[i,"type"]), eps = hs[i])
                        }
                }
                siafInts
            }
        # Matrix of terms in the epidemic component
        eTerms <- cbind(
            qSum = qSumTypes[eventData$type],
            expeta = exp(drop(mme %*% gamma)),
            siafInt = siafInt
        )
        # Return
        list(eTerms, .bdist, .influenceRegion)
    }


    ### function calculating the (upper bound) intensity of the ground process
    ### it relies on several objects for the epidemic component which are updated alongside simulation

    # t will be one of the break points in stgrid or an event time
    lambdagVec <- function (t, upper=FALSE)
    {
        ## endemic part
        hIntWKt <- hIntWK[[as.character(tBLOCK)]]

        ## epidemic part
        if (length(infectives) == 0L) ejIntWt <- numeric(0L) else {
            eTerms <- eTerms[infectives,,drop=FALSE]
            gTerm <- if (upper) {
                    rep.int(gmax, length(infectives))
                } else {
                    times <- eventMatrix[infectives,"time"]
                    types <- eventMatrix[infectives,"type"]
                    tiaf$g(t-times, tiafpars, types)
                }
            # ejIntWt only for infectives, others have 0
            ejIntWt <- apply(cbind(eTerms,gTerm), 1, prod)
            names(ejIntWt) <- infectives
        }

        c("0"=hIntWKt, ejIntWt)   # endemic component has index "0" !
    }


    ### helper function calculating the integral of lambdag from oldct to ct
    ### during simulation; it depends on the current values of the simulation

    add2Lambdag <- if (constanttiaf) {
            function () lambdagUpper * (ct-oldct)
        } else function () {
            # old endemic ground intensity * passed time
            hIntWKInt_oldct_ct <- lambdaghe[1L] * (ct-oldct)
            # integrated epidemic ground intensities of infectives (from oldct)
            ejIntWInt_oldct_ct <- if (length(infectives) == 0L) numeric(0L) else {
                eTermsProd <- apply(eTerms[infectives,,drop=FALSE], 1, prod)
                # integral of \id_{(0;eps.t]}(t-t_j) g(t-t_j \vert \kappa_j) from oldct to ct, for j in infectives
                # we can ignore the indicator because t-t_j is not >eps.t if t in [oldct;ct], because recoveries are change points
                times <- eventMatrix[infectives,"time"]
                types <- eventMatrix[infectives,"type"]
                gInt_0_ct    <- tiaf$G(ct   -times, tiafpars, types)
                gInt_0_oldct <- tiaf$G(oldct-times, tiafpars, types)
                gInt_oldct_ct <- gInt_0_ct - gInt_0_oldct
                eTermsProd * gInt_oldct_ct
            }
            sum(hIntWKInt_oldct_ct, ejIntWInt_oldct_ct)
        }



    ##################
    ### Simulation ###
    ##################

    ### Initialise values for simulation loop

    # all necessary components for an epidataCS object will be build along the simulation
    # let's start with the events of the prehistory
    tmp <- eTermsCalc(eventData, eventCoords)
    eTerms <- tmp[[1]]; rownames(eTerms) <- NULL
    bdists <- tmp[[2]]
    influenceRegions <- tmp[[3]]
    sources <- lapply(seq_len(Nout), function(x) integer(0L))

    # Transform eventData into a matrix, which is faster with rbind
    # (factors will be recreated at the end of simulation)
    # simulated events will be subsequently appended to this matrix
    eventMatrix <- if (Nout == 0L) {
        matrix(numeric(0L), nrow=0L, ncol=ncol(eventData), dimnames=list(NULL, names(eventData)))
    } else {
        sapply(eventData, as.numeric, simplify = TRUE)    # prehistory
    }
    if (Nout == 1L) eventMatrix <- t(eventMatrix)
    # we will also know about the source of infection and corresponding BLOCK in stgrid
    navec <- rep.int(NA_real_, Nout)
    eventMatrix <- cbind(eventMatrix, source = navec, lambda.h = navec,
                         lambda.e = navec, Lambdag = navec, BLOCK = navec)

    # row indices of currently infective individuals
    infectives <- seq_len(Nout)

    # maximum total number of events (including prehistory)
    maxEvents <- Nout + nEvents

    # change points of lambdag
    stgridbreaks <- blockstarts[-1,2]
    Rtimes <- structure(eventMatrix[,"time"]+eventMatrix[,"eps.t"], names=eventMatrix[,"ID"])

    # index of next event (row in eventMatrix)
    j <- Nout + 1L

    # allocation of large objects for faster filling-in of new events
    allocated <- Nout
    ncolEventMatrix <- ncol(eventMatrix)
    newAllocation <- expression({
        eventMatrix <- rbind(eventMatrix, matrix(NA_real_, nrow = .allocate, ncol = ncolEventMatrix))
        eventCoords <- rbind(eventCoords, matrix(NA_real_, nrow = .allocate, ncol = 2L))
        eTerms <- rbind(eTerms, matrix(NA_real_, nrow = .allocate, ncol = 3L))
        bdists <- c(bdists, rep.int(NA_real_,.allocate))
        influenceRegions <- c(influenceRegions, vector(.allocate, mode="list"))
        sources <- c(sources, vector(.allocate, mode="list"))
        allocated <- allocated + .allocate
    })

    # current time point
    ct <- t0

    # current value of the cumulative intensity function of the ground process
    Lambdag <- 0

    # last point rejected?
    pointRejected <- FALSE

    # did we have numerical problems simulating from Exp(lambdagUpper) in the current loop?
    hadNumericalProblemsInf <- hadNumericalProblems0 <- FALSE

    # index of the current loop
    loopCounter <- 0L


    ### Let's Rock 'n' Roll

    if (trace > 0L) {
        cat("\nSimulation path (starting from t=", t0, "):\n---\n", sep="")
    } else {
        cat("\nSimulating (starting from t=", t0, ") ...\n", sep="")
    }

    while(j <= maxEvents && ct < T)
    {
        loopCounter <- loopCounter + 1L
        if (trace > 0L && loopCounter %% trace == 0L) {
            cat(loopCounter, "@t =", ct, ":\t#simulated events =", j-1L-Nout,
            "\t#currently infective =", length(infectives),
            if (hase && !constanttiaf) paste("\tlast rejected?", pointRejected), "\n")
            flush.console()   # affects Windows only
        }

        # check if we need to allocate larger matrices
        if (j > allocated) {
            eval(newAllocation)
        }

        if (!pointRejected)  # what we have to do in the usual case
        {
            # we need the time block of stgrid corresponding to the new covariates,
            # i.e. search BLOCK such that t in [start; stop)
            tBLOCK <- blockstarts[findInterval(ct, blockstarts[,2]), 1]

            # Compute new infection intensity (upper bound)
            lambdaghe <- lambdagVec(ct, upper=TRUE)
            lambdagUpper <- sum(lambdaghe)

            # Determine time of next external change point
            changePoints <- c(nextblock = if (length(stgridbreaks) > 0L) stgridbreaks[1L],
                              Rtimes)
            nextChangePoint <- if (length(changePoints) > 0L) {
                   changePoints[which.min(changePoints)]    # don't use min() because need names
                } else Inf
        }
        pointRejected <- FALSE

        ## Simulate waiting time for the subsequent infection
        Delta <- tryCatch(rexp(1, rate = lambdagUpper),
            warning = function (w) {
                if (lambdagUpper < 1) {   # rate was to small for rexp
                    assign("hadNumericalProblems0", TRUE, inherits = TRUE)
                    if (nextChangePoint == Inf) NULL else Inf
                } else {   # rate was to big for rexp
                    0   # since R-2.7.0 rexp(1, Inf) returns 0 with no warning!
                }
            })

        # Stop if lambdaStarMax too small AND no more changes in rate
        if (is.null(Delta)) break
        # Stop if lambdaStarMax too big meaning Delta == 0 (=> concurrent events)
        if (Delta == 0) {
            hadNumericalProblemsInf <- TRUE
            break
        }
        # Stop at all costs if end of simulation time [t0; T) has been reached
        if (isTRUE(min(ct+Delta, nextChangePoint) >= T)) break
        # ">=" because we don't want an event at "end"

        oldct <- ct
        if (ct + Delta > nextChangePoint) {
        ## Simulated time point is beyond the next time of intensity change (removal or endemic covariates)
            ct <- unname(nextChangePoint)
            # update cumulative intensity of the ground processes up to time ct,
            # i.e. add integral of lambdag from oldct to ct
            Lambdag <- Lambdag + add2Lambdag()
            # is this change point due to next time block in stgrid?
            if (names(nextChangePoint) == "nextblock") {
                stgridbreaks <- stgridbreaks[-1]
            } else { # i.e. change point due to recovery
                recoverer <- names(nextChangePoint)
                # update set of infectives
                infectives <- setdiff(infectives, recoverer)
                # remove recovery time from Rtimes
                .Rtimesidx <- match(recoverer, names(Rtimes))
                Rtimes <- Rtimes[-.Rtimesidx]
            }
        } else {
        ## Simulated time point lies within the thinning period
            ct <- ct + Delta
            # rejection sampling if non-constant temporal interaction kernel g
            if (hase && !constanttiaf) {
                # Calculate actual ground intensity for rejection probability at new ct
                lambdaghe <- lambdagVec(ct, upper=FALSE)
                lambdag <- sum(lambdaghe)
                # rejection sampling step
                if (lambdag/lambdagUpper < runif(1)) {
                    pointRejected <- TRUE
                    next
                }
            }
            # At this point, we have an actual event!

            # update cumulative intensity of the ground processes up to time ct,
            # i.e. add integral of lambdag from oldct to ct
            Lambdag <- Lambdag + add2Lambdag()
            # note that lambdaghe[1L] did not change by the above update in case of !constanttiaf,
            # which is expected by add2Lambdag (which requires the value of lambdag.h(oldct))

            # Where did the event come from: imported case or infection?
            .eventSource <- as.integer(sample(names(lambdaghe), 1L, prob=lambdaghe))

            # We now sample type and location
            if (.eventSource == 0L) {   # i.e. endemic source of infection
                .eventType <- sample(typeNames, 1L, prob=if (nbeta0 > 1L) exp(beta0))
                stgrididx <- which(gridBlocks == tBLOCK)
                .eventTile <- sample(stgrid$tile[stgrididx], 1L, prob=dsexpeta[stgrididx])  # this is a factor
                # in spsample it is not guaranteed that the sample will consists of exactly n=1 point
                while(inherits(
                    eventLocationSP <- try(spsample(tiles[as.character(.eventTile),], n=1L, type="random"), silent = TRUE)
                , "try-error")) {}  # if no point is sampled (very unlikely though) we would get an error
                .eventLocation <- coordinates(eventLocationSP)[1L,,drop=FALSE]
            } else {    # i.e. source is one of the currently infective individuals
                sourceType <- eventMatrix[.eventSource,"type"]
                sourceCoords <- eventCoords[.eventSource,,drop=FALSE]
                sourceIR <- influenceRegions[[.eventSource]]
                .eventType <- sample(typeNames[qmatrix[sourceType,]], 1L)
                .eventTypeCode <- match(.eventType, typeNames)
                eventLocationIR <- if (constantsiaf) {
                        as.matrix(spatstat::coords(spatstat::runifpoint(1L, win=sourceIR, giveup=1000)))
                    } else {
                        eventInsideIR <- FALSE
                        ntries <- 0L
                        while(!eventInsideIR) {
                            if (ntries >= 1000) { # TODO: throw warning, skip to next event?
                                cat("rejection sampling of event location seems difficult...\n")
                                browser()
                            }
                            ntries <- ntries + 1L
                            eventLocationIR <- siaf$simulate(1L, siafpars, .eventTypeCode)
                            eventInsideIR <- spatstat::inside.owin(eventLocationIR[,1],
                                  eventLocationIR[,2], sourceIR)
                        }
                        eventLocationIR
                    }
                .eventLocation <- sourceCoords + eventLocationIR
                whichTile <- overlay(SpatialPoints(.eventLocation, proj4string=CRS(proj4string(W))), tiles)
                .eventTile <- row.names(tiles)[whichTile]
                .eventTile <- factor(.eventTile, levels=tileLevels)
            }
            .eventType <- factor(.eventType, levels=typeNames)

            # sample marks at this time and location
            .eventMarks <- rmarks(ct, .eventLocation)

            # gather event information
            .eventData <- data.frame(ID=j, time=ct, tile=.eventTile, type=.eventType,
                .eventMarks, check.rows = FALSE, check.names = FALSE)

            # determine potential sources of infection (for epidataCS and lambda)
            .sources <- infectives[eventMatrix[infectives,"type"] %in% which(qmatrix[,.eventType])]
            if (length(.sources) > 0L) {
                .sdiffs <- .eventLocation[rep.int(1L,length(.sources)),,drop=FALSE] - eventCoords[.sources,,drop=FALSE]
                .sources <- .sources[sqrt(rowSums(.sdiffs^2)) <= eventMatrix[.sources,"eps.s"]]
            }

            # calculate actual intensity at this time, location and type
            .mmhEvent <- buildmmh(.eventData)
            .etaEvent <- .mmhEvent %*% beta
            .offsetEvent <- attr(.mmhEvent, "offset")
            if (!is.null(.offsetEvent)) .etaEvent <- .offsetEvent + .etaEvent
            .lambdah <- exp(.etaEvent)
            .lambdae <- if (length(.sources) > 0L) {
                .sdiffs <- .eventLocation[rep.int(1L,length(.sources)),,drop=FALSE] - eventCoords[.sources,,drop=FALSE]
                .fSources <- siaf$f(.sdiffs, siafpars, eventMatrix[.sources,"type"])
                .gSources <- tiaf$g(ct - eventMatrix[.sources,"time"], tiafpars, eventMatrix[.sources,"type"])
                sum(eTerms[.sources,"expeta"] * .fSources * .gSources)
            } else 0

            # calculate terms of the epidemic component e_j(t,s) of the new infective
            tmp <- eTermsCalc(.eventData, .eventLocation)

            # Update objects
            eventMatrix[j,] <- c(j, ct, as.numeric(.eventTile), as.numeric(.eventType),
                                 sapply(.eventMarks, as.numeric), .eventSource,
                                 .lambdah, .lambdae, Lambdag, tBLOCK)
            eventCoords[j,] <- .eventLocation
            eTerms[j,] <- tmp[[1]]
            bdists[j] <- tmp[[2]]
            influenceRegions[[j]] <- tmp[[3]][[1]]
            sources[[j]] <- .sources

            # Update set of infectives and recovery times
            infectives <- c(infectives, j)
            Rtimes <- c(Rtimes, structure(ct + .eventMarks[["eps.t"]], names = j))

            # Increment next event iterator
            j <- j + 1L
        }
    }
    if (trace > 0L) cat("---\n")

    
    ### update T if simulation ended because of limited 'nEvents'

    if (j > maxEvents) {
        T <- ct
        # clip stgrid to effective time range of simulation
        stgrid <- subset(stgrid, start <= T)
    }
    cat("Simulation has ended @t =", T, "with", j-1L-Nout, "simulated events.\n")



    ##############
    ### Return ###
    ##############


    ### Throw warnings in case of numerical difficulties

    if (hadNumericalProblemsInf) {
        warning("simulation ended due to an infinite overall infection rate")
    }
    if (hadNumericalProblems0) {
        warning("occasionally, the overall infection rate was numerically equal to 0")
    }


    ### transform eventMatrix back into a data.frame with original factor variables

    if (trace > 0L) {
      cat("\nConverting simulated events into an object of class \"epidataCS\"...\n")
    }
    preEventData <- eventData

    # drop unused entries (due to large pre-allocation) from objects
    seqAlongEvents <- seq_len(j-1L)
    eventData <- as.data.frame(eventMatrix[seqAlongEvents,,drop=FALSE])

    # rebuild factor variables
    for (idx in which(sapply(preEventData, is.factor))) {
        origlevels <- levels(preEventData[[idx]])
        eventData[[idx]] <- factor(eventData[[idx]], levels=seq_along(origlevels), labels=origlevels)
    }

    # transform integer columns to integer
    eventData[c("ID","source","BLOCK")] <- lapply(eventData[c("ID","source","BLOCK")], as.integer)


    ### Append additional columns for an epidataCS object

    # add endemic covariates at events
    stgrididx <- apply(eventData[c("BLOCK","tile")], 1, function (x) {
        ret <- with(stgrid, which(BLOCK==as.integer(x[1L]) & tile==x[2L]))
        if (length(ret) == 0L) NA_integer_ else ret
        #<- events of the prehistory have missing BLOCKs, thus return NA
    })
    stgridIgnoreCols <- match(c("BLOCK",setdiff(obligColsNames_stgrid, "start")), names(stgrid))
    eventData <- cbind(eventData, stgrid[stgrididx, -stgridIgnoreCols])
    rownames(eventData) <- eventData$ID

    # add hidden columns
    eventData$.obsInfLength <- with(eventData, pmin(T-time, eps.t))
    eventData$.sources <- sources[seqAlongEvents]
    eventData$.bdist <- bdists[seqAlongEvents]
    eventData$.influenceRegion <- influenceRegions[seqAlongEvents]


    ### Construct "epidataCS" object

    events <- SpatialPointsDataFrame(
        coords = eventCoords[seqAlongEvents,,drop=FALSE],
        data = eventData,
        proj4string = CRS(proj4string(W)),
        match.ID = FALSE
    )

    if (.onlyEvents) {
        cat("Done.\n")
        attr(events, "timeRange") <- c(t0, T)
        attr(events, "runtime") <- proc.time()[[3]] - ptm
        return(events)
    }

    epi <- list(events=events, stgrid=stgrid, W=W, qmatrix=qmatrix)


    ### Return object of class "simEpidataCS"

    cat("Done.\n")
    # append configuration of the model
    epi$timeRange <- c(t0, T)
    epi$formula <- list(
                   endemic = if (typeSpecificEndemicIntercept) {
                       update(formula(endemic), ~ (1|type) + .)   # re-add to the formula
                   } else formula(endemic),
                   epidemic = formula(epidemic),
                   siaf = siaf, tiaf = tiaf
                   )
    epi$coefficients <- list(beta0=beta0, beta=beta, gamma=gamma,
                             siafpars=siafpars, tiafpars=tiafpars)
    epi$npars <- c(nbeta0=nbeta0, p=p, q=q, nsiafpars=nsiafpars, ntiafpars=ntiafpars)
    epi$call <- cl
    epi$runtime <- proc.time()[[3]] - ptm
    class(epi) <- c("simEpidataCS", "epidataCS", "list")
    return(epi)
}






################################################################################
# A 'simulate' method for objects of class "twinstim".
################################################################################

### TODO: actually stgrid's of simulations might have different time ranges when nEvents is active -> atm, simplify ignores this

simulate.twinstim <- function (object, nsim = 1, seed = NULL, data, tiles,
    rmarks = NULL, t0 = NULL, T = NULL, nEvents = 1e5, nCub,
    W = NULL, trace = FALSE, nCircle2Poly = 32, gmax = NULL, .allocate = 500, simplify = TRUE, ...)
{
    ptm <- proc.time()[[3]]
    cl <- match.call()


    ### Determine seed (copied from stats:::simulate.lm)

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }


    ### Few checks

    stopifnot(inherits(object, "twinstim"), inherits(data, "epidataCS"))
    stopifnot(isScalar(nsim), nsim > 0)
    nsim <- as.integer(nsim)
    if (is.null(t0)) t0 <- object$timeRange[1]
    if (is.null(T))  T  <- object$timeRange[2]


    ### Retrieve arguments for simulation

    endemic  <- formula(object)$endemic
    epidemic <- formula(object)$epidemic
    # we don't need any reference to the original twinstim evaluation environment
    environment(endemic) <- environment(epidemic) <- globalenv()
    siaf     <- formula(object)$siaf
    tiaf     <- formula(object)$tiaf
    if (is.null(rmarks)) {
        observedMarks <- subset(marks(data), subset = time > t0 & time <= T)
        observedMarks <- observedMarks[match("eps.t", names(observedMarks)):(ncol(observedMarks)-2L)]
        rmarks <- function (t, s) {
            as.data.frame(lapply(observedMarks, function (x) sample(na.omit(x), size=1L)), optional=TRUE)
        }
    }
    theta <- coef(object)
    beta0 <- beta <- gamma <- siafpars <- tiafpars <- NULL
    with(as.list(object$npars), {
        beta0    <<- theta[seq_len(nbeta0)]
        beta     <<- theta[nbeta0+seq_len(p)]
        gamma    <<- theta[nbeta0+p+seq_len(q)]
        siafpars <<- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <<- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]
    })


    ### Run the simulation(s)

    # First simulation
    if (nsim > 1L) {
        cat("\nTime at beginning of simulation:", as.character(Sys.time()), "\n")
        cat("Simulation 1 /", nsim, "...\n")
        cat("-------------------------------------------------------------------------------\n")
    }
    res <- simEpidataCS(endemic, epidemic, siaf, tiaf, object$qmatrix, rmarks,
                        data$events, data$stgrid, tiles, beta0, beta, gamma, siafpars,
                        tiafpars, t0, T, nEvents, nCub, W, trace, nCircle2Poly,
                        gmax, .allocate, .skipChecks = TRUE)
    if (nsim > 1L) {
        cat("\n-------------------------------------------------------------------------------\n")
        cat("Runtime of first simulation:", res$runtime, "seconds\n")
        cat("Estimated finishing time:", as.character(Sys.time() + (nsim-1) * res$runtime), "\n\n")
        # set up list of simulations
        res <- if (simplify) {
                with(res, list(
                    eventsList=c(structure(events, timeRange = timeRange, runtime = runtime),
                                 vector(nsim-1L, mode="list")),
                    stgrid=stgrid, W=W, qmatrix=qmatrix, formula=formula,
                    coefficients=coefficients, npars=npars, call=call
                ))
            } else {
                c(list(res), vector(nsim-1L, mode="list"))
            }
        # force garbage collection
        capture.output(gc())
        # run the remaining simulations
        for (i in 2:nsim) {
            cat("Simulation", sprintf(paste0("%",nchar(nsim),"i"), i), "/", nsim, "...")
            capture.output(
            resi <- simEpidataCS(endemic, epidemic, siaf, tiaf, object$qmatrix, rmarks,
                data$events, data$stgrid, tiles, beta0, beta, gamma, siafpars,
                tiafpars, t0, T, nEvents, nCub, W, trace, nCircle2Poly,
                gmax, .allocate, .skipChecks = TRUE, .onlyEvents = simplify)
            )
            cat("\tsimulated", if (simplify) sum(!is.na(resi$source)) else sum(!is.na(resi$events$source)),
                "events up to time", if (simplify) attr(resi,"timeRange")[2] else resi$timeRange[2], "\n")
            if (simplify) res$eventsList[[i]] <- resi else res[[i]] <- resi
        }
        cat("\nDone (", as.character(Sys.time()), ").\n", sep="")
    }
    attr(res, "call") <- cl
    attr(res, "seed") <- RNGstate
    attr(res, "simplified") <- simplify
    attr(res, "runtime") <- proc.time()[[3]] - ptm
    class(res) <- if (nsim == 1L) {
            c("simEpidataCS", "epidataCS", "list")
        } else c("simEpidataCSlist", "list")
    res
}



### print method for lists of simulated epidemics

print.simEpidataCSlist <- function (x, ...)
{
    cat("\nCall:\n")
    print.default(attr(x, "call"))
    simplified <- attr(x, "simplified")
    nsim <- if (simplified) length(x$eventsList) else length(x)
    cat("\n")
    cat(if (simplified) "Simplified list" else "List", "of", nsim,
        "simulated epidemics of class \"simEpidataCS\" (not printed)\n\n")
    invisible(x)
}
