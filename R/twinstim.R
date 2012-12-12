################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Maximum Likelihood inference for the two-component spatio-temporal intensity
### model described in Meyer et al (2012), DOI: 10.1111/j.1541-0420.2011.01684.x
###
### Copyright (C) 2012 Sebastian Meyer
### $Revision: 454 $
### $Date: 2012-11-15 15:43:56 +0100 (Do, 15. Nov 2012) $
################################################################################


twinstim <- function (endemic, epidemic, siaf, tiaf, qmatrix = data$qmatrix,
    data, subset, t0 = data$stgrid$start[1], T = tail(data$stgrid$stop,1),
    na.action = na.fail, partial = FALSE,
    control.siaf = list(F=list(), Deriv=list()), optim.args, finetune = FALSE,
    model = FALSE, cumCIF = TRUE, cumCIF.pb = TRUE)
{

    ####################
    ### Preparations ###
    ####################

    ptm <- proc.time()[[3]]
    cl <- match.call()
    partial <- as.logical(partial)
    finetune <- if (partial) FALSE else as.logical(finetune)

    ## # Collect polyCub.midpoint warnings (and apply unique on them at the end)
    ## .POLYCUB.WARNINGS <- NULL
    ## .polyCub <- function (...) {
    ##     int <- withCallingHandlers(
    ##            polyCub.midpoint(...),
    ##            warning = function (w) {
    ##                .POLYCUB.WARNINGS <<- c(.POLYCUB.WARNINGS, list(w))
    ##                invokeRestart("muffleWarning")
    ##            })
    ## }
    
    # Clean the model environment when exiting the function
    on.exit(suppressWarnings(rm(cl, cumCIF, cumCIF.pb, data, doHessian,
        eventDists, eventsData, finetune, neghess, fisherinfo, fit, fixed,
        functions, globalEndemicIntercept, inmfe, initpars, ll, negll, loglik,
        mfe, mfhEvents, mfhGrid, model, my.na.action, na.action, namesOptimUser,
        namesOptimArgs, nlminbControl, nlminbRes, nlmObjective, nlmControl,
        nlmRes, nmRes, optim.args, optimArgs, control.siaf,
        optimMethod, optimRes, optimRes1, optimValid, partial,
        partialloglik, ptm, qmatrix, res, negsc, score, subset, tmpexpr,
        typeSpecificEndemicIntercept, useScore, whichfixed, 
        inherits = FALSE)))


    ### Verify that 'data' inherits from "epidataCS"

    if (!inherits(data, "epidataCS")) {
        stop("'data' must inherit from class \"epidataCS\"")
    }


    ### Check time range

    if (!isScalar(t0) || !isScalar(T)) {
        stop("endpoints 't0' and 'T' must be single numbers")
    }
    if (T <= t0) {
        stop("'T' must be greater than 't0'")
    }
    if (!t0 %in% data$stgrid$start) {
        justBeforet0 <- match(TRUE, data$stgrid$start > t0) - 1L
        # if 't0' is beyond the time range covered by 'data$stgrid'
        if (is.na(justBeforet0)) justBeforet0 <- length(data$stgrid$start)   # t0 was too big
        if (justBeforet0 == 0L) justBeforet0 <- 1L   # t0 was too small
        t0 <- data$stgrid$start[justBeforet0]
        message("replaced 't0' by the value ", t0,
                " (must be a 'start' time of 'data$stgrid')")
    }
    if (!T %in% data$stgrid$stop) {
        justAfterT <- match(TRUE, data$stgrid$stop > T)
        # if 'T' is beyond the time range covered by 'data$stgrid'
        if (is.na(justAfterT)) justAfterT <- length(data$stgrid$stop)   # T was too big
        T <- data$stgrid$stop[justAfterT]
        message("replaced 'T' by the value ", T,
                " (must be a 'stop' time of 'data$stgrid')")
    }


    ### Subset events

    eventsData <- if (missing(subset)) data$events@data else {
        do.call("subset.data.frame", args = list(
            x = quote(data$events@data), subset = cl$subset, drop = FALSE
        ))
    }






    #############################################################
    ### Build up a model.frame for both components separately ###
    #############################################################



    ##########################
    ### epidemic component ###
    ##########################


    ### Parse epidemic formula

    if (missing(epidemic)) {
        epidemic <- ~ 0
    } else {
        environment(epidemic) <- environment()
        # so that t0 and T are found in the subset expressions below
    }
    epidemic <- terms(epidemic, data = eventsData, keep.order = TRUE)
    if (!is.null(attr(epidemic, "offset"))) {
        warning("offsets are not implemented for the 'epidemic' component")
    }


    ### Generate model frame

    # na.action mod such that for simulated epidataCS, where events of the
    # prehistory have missing 'BLOCK' indexes, those NA's do not matter.
    # ok because actually, 'eventBlocks' are only used in the partial likelihood
    # and there only eventBlocks[includes] is used (i.e. no prehistory events)
    my.na.action <- function (object, ...) {
        prehistevents <- object[object[["(time)"]] <= t0, "(ID)"]
        if (length(prehistevents) == 0L) return(na.action(object, ...))
        origprehistblocks <- object[match(prehistevents,object[["(ID)"]]), "(BLOCK)"]
        object[object[["(ID)"]] %in% prehistevents, "(BLOCK)"] <- 0L
        xx <- na.action(object, ...)
        xx[match(prehistevents,xx[["(ID)"]],nomatch=0L), "(BLOCK)"] <-
            origprehistblocks[prehistevents %in% xx[["(ID)"]]]
        xx
    }

    ID <- tile <- type <- BLOCK <- .obsInfLength <- .bdist <-
        "just cheating on codetools::checkUsage"
    mfe <- model.frame(epidemic, data = eventsData,
                       subset = time + eps.t > t0 & time <= T,
# here we can have some additional rows (individuals) compared to mfhEvents, which is established below!
# Namely those with time in (t0-eps.t; t0], i.e. still infective individuals, which are part of the prehistory of the process
                       na.action = my.na.action,
# since R 2.10.0 patched also works with epidemic = ~1 and na.action=na.fail (see PR#14066)
                       drop.unused.levels = FALSE,
                       ID = ID, time = time, tile = tile, type = type,
                       eps.t = eps.t, eps.s = eps.s, BLOCK = BLOCK,
                       obsInfLength = .obsInfLength, bdist = .bdist)
    rm(ID, tile, type, BLOCK, .obsInfLength, .bdist)

    ### Extract essential information from model frame

    # inmfe=rowindex(data$events@data) is necessary for subsetting
    # influenceRegion (list object not compatible with model.frame) and
    # coordinates
    # CAVE: ID is not necessarily identical to the rowindex of data$events,
    #       e.g., if working with subsetted epidataCS
    inmfe <- which(data$events@data$ID %in% mfe[["(ID)"]])
    N <- length(inmfe)   # mfe also contains events of the prehistory
    eventTimes <- mfe[["(time)"]] # I don't use model.extract since it returns named vectors
    # Indicate events after t0, which are actually part of the process
    # (events in (-Inf;t0] only contribute in sum over infected individuals)
    includes <- which(eventTimes > t0)   # this indexes mfe!
    Nin <- length(includes)
    if (Nin == 0L) {
        stop("none of the ", nrow(data$events@data), " supplied ",
             "events is in the model (check 'subset', 't0' and 'T')")
    }
    eventBlocks <- mfe[["(BLOCK)"]]   # only necessary for partial log-likelihood
    eventTypes <- factor(mfe[["(type)"]])   # drop unused levels
    typeNames <- levels(eventTypes)
    nTypes <- length(typeNames)
    if (nTypes > 1L) cat("marked point pattern of", nTypes, "types\n")
    qmatrix <- checkQ(qmatrix, typeNames)
    # we only need the integer codes for the calculations
    eventTypes <- as.integer(eventTypes)


    ### Generate model matrix

    mme <- model.matrix(epidemic, mfe)
    q <- ncol(mme)
    hase <- q > 0L


    ### Extract further model components (only if q > 0)

    if (hase) {
        eps.t <- mfe[["(eps.t)"]]
        removalTimes <- eventTimes + eps.t
        eps.s <- mfe[["(eps.s)"]]
        bdist <- mfe[["(bdist)"]]
        gIntUpper <- mfe[["(obsInfLength)"]]
        gIntLower <- pmax(0, t0-eventTimes)
        eventCoords <- coordinates(data$events)[inmfe,,drop=FALSE]
        eventDists <- as.matrix(dist(eventCoords, method = "euclidean"))
        influenceRegion <- data$events@data$.influenceRegion[inmfe]
        iRareas <- sapply(influenceRegion, attr, "area")
        # determine possible event sources (need to re-do this because
        # subsetting has crashed old row indexes from the epidataCS object)
        # actually only eventSources of includes are needed
        eventSources <- lapply(1:N, function (i) {
            determineSources(i, eventTimes, removalTimes, eventDists[i,], eps.s, eventTypes, qmatrix)
        })
        # calculate sum_{k=1}^K q_{kappa_j,k} for all j = 1:N
        qSum <- unname(rowSums(qmatrix)[eventTypes])   # N-vector
    } else message("no epidemic component in model")



    #########################
    ### endemic component ###
    #########################


    ### Parse endemic formula

    if (missing(endemic)) {
        endemic <- ~ 0
    } else {
        environment(endemic) <- environment()
        # so that t0 and T are found in the subset expressions below
    }
    endemic <- terms(endemic, data = data$stgrid, keep.order = TRUE)

    ## check for type-specific endemic intercept and remove it from the formula
    ## (will be handled separately)
    typeSpecificEndemicIntercept <- "1 | type" %in% attr(endemic, "term.labels")
    if (typeSpecificEndemicIntercept) {
        endemic <- update(endemic, ~ . - (1|type)) # this drops the terms attributes
        endemic <- terms(endemic, data = data$stgrid, keep.order = TRUE)
    }

    globalEndemicIntercept <- if (typeSpecificEndemicIntercept) {
            attr(endemic, "intercept") <- 1L   # we need this to ensure that we have correct contrasts
            FALSE
        } else attr(endemic, "intercept") == 1L

    nbeta0 <- globalEndemicIntercept + typeSpecificEndemicIntercept * nTypes


    ### Generate endemic model frame and model matrix on event data

    ID <- "just cheating on codetools::checkUsage"
    mfhEvents <- model.frame(endemic, data = eventsData,
                             subset = time>t0 & time<=T & ID %in% mfe[["(ID)"]],
                             na.action = na.fail,
                             # since R 2.10.0 patched also works with
                             # endemic = ~1 (see PR#14066)
                             drop.unused.levels = FALSE)
    rm(ID)
    mmhEvents <- model.matrix(endemic, mfhEvents)
    # exclude intercept from endemic model matrix below, will be treated separately
    if (nbeta0 > 0) mmhEvents <- mmhEvents[,-1,drop=FALSE]
    #stopifnot(nrow(mmhEvents) == Nin)
    p <- ncol(mmhEvents)
    hash <- (nbeta0+p) > 0L


    ### Generate model frame and model matrix on grid data (only if p > 0)

    if (hash) {
        offsetEvents <- model.offset(mfhEvents)
        BLOCK <- tile <- area <- "just cheating on codetools::checkUsage"
        mfhGrid <- model.frame(endemic, data = data$stgrid,
                               subset = start >= t0 & stop <= T,
                               na.action = na.fail,
                               # since R 2.10.0 patched also works with
                               # endemic = ~1 (see PR#14066)
                               drop.unused.levels = FALSE,
                               BLOCK=BLOCK, tile=tile, dt=stop-start, ds=area)
                               # 'tile' is redundant here for fitting but useful
                               # for debugging & necessary for intensityplots
        rm(BLOCK, tile, area)
        gridBlocks <- mfhGrid[["(BLOCK)"]]
        histIntervals <- unique(data$stgrid[c("BLOCK", "start", "stop")]) # sorted
        rownames(histIntervals) <- NULL
        histIntervals <- histIntervals[histIntervals$start >= t0 &
                                       histIntervals$stop <= T,]
        gridTiles <- mfhGrid[["(tile)"]] # only needed for intensityplot
        mmhGrid <- model.matrix(endemic, mfhGrid)

        # exclude intercept from endemic model matrix below, will be treated separately
        if (nbeta0 > 0) mmhGrid <- mmhGrid[,-1,drop=FALSE]
        # Extract endemic model components
        offsetGrid <- model.offset(mfhGrid)
        dt <- mfhGrid[["(dt)"]]
        ds <- mfhGrid[["(ds)"]]
    } else message("no endemic component in model")


    ### Check that there is at least one parameter

    if (!hash && !hase) {
        stop("nothing to do: neither endemic nor epidemic parts were specified")
    }






    #############################
    ### Interaction functions ###
    #############################

    if (hase) {

        ## Check interaction functions
        siaf <- do.call(".parseiaf", args = alist(siaf))
        constantsiaf <- attr(siaf, "constant")
        nsiafpars <- siaf$npars

        tiaf <- do.call(".parseiaf", args = alist(tiaf))
        constanttiaf <- attr(tiaf, "constant")
        ntiafpars <- tiaf$npars

        ## Check control.siaf
        if (constantsiaf) control.siaf <- NULL else {
            stopifnot(is.null(control.siaf) || is.list(control.siaf))
            if (is.list(control.siaf)) stopifnot(sapply(control.siaf, is.list))
        }
        
        ## Define function that integrates the 'tiaf' function
        .tiafInt <- .tiafIntFUN()

        ## Define function that integrates the two-dimensional 'siaf' function
        ## over the influence regions of the events
        .siafInt <- .siafIntFUN(siaf = siaf, noCircularIR = all(eps.s > bdist))
        .siafInt.args <- c(alist(siafpars), control.siaf$F)

        ## Memoisation of .siafInt
        ..siafInt <- if (requireNamespace("memoise")) {
            memoise::memoise(.siafInt)
            ## => speed-up optimization since 'nlminb' evaluates the loglik and
            ## score for the same set of parameters at the end of each iteration
        } else {
            cat("Continuing without memoisation of 'siaf$f' cubature ...\n")
            ## However, trivial caching is used manually in case of fixed
            ## "siafpars" and in LambdagEvents()
            .siafInt
        }

    } else {
        
        if (!missing(siaf) && !is.null(siaf))
            warning("'siaf' can only be modelled in conjunction with an 'epidemic' process")
        if (!missing(tiaf) && !is.null(tiaf))
            warning("'tiaf' can only be modelled in conjunction with an 'epidemic' process")
        siaf <- tiaf <- NULL
        nsiafpars <- ntiafpars <- 0L
        control.siaf <- NULL

    }

    hassiafpars <- nsiafpars > 0L
    hastiafpars <- ntiafpars > 0L

    ## Can we calculate the score function?
    useScore <- if (partial) FALSE else if (hase) {
        (!hassiafpars | !is.null(siaf$deriv)) &
        (!hastiafpars | (!is.null(tiaf$deriv)) & !is.null(tiaf$Deriv))
    } else TRUE

    ## Define function that applies siaf$Deriv on all events (integrate the
    ## two-dimensional siaf$deriv function)
    if (useScore && hassiafpars) {
        .siafDeriv <- function (siafpars, ...) {
            derivInti <- function (i)
                siaf$Deriv(influenceRegion[[i]], siaf$deriv, siafpars,
                           eventTypes[i], ...)
            derivInt <- sapply(1:N, derivInti)
            #<- N-vector or nsiafpars x N matrix => transform to N x nsiafpars
            if (is.matrix(derivInt)) t(derivInt) else as.matrix(derivInt)
        }
        .siafDeriv.args <- c(alist(siafpars), control.siaf$Deriv)
    }




    ############################################################################
    ### Log-likelihood function, score function, expected Fisher information ###
    ############################################################################


    ### Total number of parameters (= length of 'theta')

    npars <- nbeta0 + p + q + nsiafpars + ntiafpars

    # REMINDER:
    #  theta - parameter vector c(beta, gamma, siafpars, tiafpars), where
    #    beta    - parameters of the endemic component exp(offset + eta_h(t,s))
    #    gamma   - parameters of the epidemic term exp(eta_e(t,s))
    #    siafpars- parameters of the epidemic spatial interaction function
    #    tiafpars- parameters of the epidemic temporal interaction function
    #  mmh[Events/Grid] - model matrix related to beta, i.e the endemic component,
    #                     either for events only or for the whole spatio-temporal grid
    #  offset[Events/Grid] - offset vector related to the endemic component (can be NULL),
    #                        either for events only or for the whole spatio-temporal grid
    #  dt, ds - columns of the spatio-temporal grid (dt = stop-start, ds = area)
    #  mme - model matrix related to gamma in the epidemic component
    #  siaf, tiaf - spatial/temporal interaction function (NULL, list or numeric)
    #  eventTimes, eventCoords, eventSources, gIntLower, gIntUpper, influenceRegion -
    #     columns of the events data frame


    if (hash)
    {
        ### Calculates the endemic component (for i in includes)
        ### h(t_i,s_i,kappa_i) = exp(offset_i + beta_{0,kappa_i} + eta_h(t_i,s_i))

        .hEvents <- function (beta0, beta) {}
        body(.hEvents) <- as.call(c(as.name("{"),
            if (p > 0L) {
                expression(eta <- drop(mmhEvents %*% beta))
            } else {
                expression(eta <- numeric(Nin))
            },
            if (nbeta0 == 1L) {
                expression(eta <- beta0 + eta)   # global intercept
            } else if (nbeta0 > 1L) {
                expression(eta <- beta0[eventTypes] + eta)   # type-specific intercept
            },
            if (!is.null(offsetEvents)) expression(eta <- offsetEvents + eta),
            expression(exp(eta))        # Nin-vector
        ))


        ### Integral of the endemic component over [0;uppert] x W

        .hIntTW <- function (beta, score = matrix(1,nrow(mmhGrid),1L), uppert = NULL) {}
        body(.hIntTW) <- as.call(c(as.name("{"),
            expression(
                subtimeidx <- if (!is.null(uppert)) { # && isScalar(uppert) && t0 <= uppert && uppert < T
                    if (uppert == t0) return(0)       # actually never happens
                                        # since uppert %in% eventTimes[includes] > t0
                    idx <- match(TRUE, histIntervals$stop >= uppert)
                    firstBlockBeyondUpper <- histIntervals$BLOCK[idx]
                    newdt <- uppert - histIntervals$start[idx]
                    dt[gridBlocks == firstBlockBeyondUpper] <- newdt
                    gridBlocks <= firstBlockBeyondUpper
                } else NULL
            ),
            if (p > 0L) {
                expression(eta <- drop(mmhGrid %*% beta))
            } else {
                expression(eta <- numeric(nrow(mmhGrid)))
            },
            if (!is.null(offsetGrid)) expression(eta <- offsetGrid + eta),
            expression(sumterms <- score * (exp(eta)*ds*dt)),
            expression(if (is.null(subtimeidx)) colSums(sumterms) else colSums(sumterms[subtimeidx,,drop=FALSE]))
        ))
    }

    if (hase)
    {
        ### Calculates the epidemic component for all events

        .eEvents <- function (gammapred, siafpars, tiafpars,
            ncolsRes = 1L, score = matrix(1,N,ncolsRes), f = siaf$f, g = tiaf$g)
            # second line arguments are for score functions with defaults for loglik
        {
            e <- matrix(0, N, ncolsRes)
            for (i in includes) {
                sources <- eventSources[[i]]
                nsources <- length(sources)
                e[i,] <- if (nsources == 0L) 0 else {
                    scoresources <- score[sources,,drop=FALSE]
                    predsources <- gammapred[sources]
                    repi <- rep.int(i, nsources)
                    sdiff <- eventCoords[repi,,drop=FALSE] - eventCoords[sources,,drop=FALSE]
                    fsources <- f(sdiff, siafpars, eventTypes[sources])
                    tdiff <- eventTimes[repi] - eventTimes[sources]
                    gsources <- g(tdiff, tiafpars, eventTypes[sources])
        # if(length(predsources) != NROW(fsources) || NROW(fsources) != NROW(gsources)) browser()
                    colSums(scoresources * predsources * fsources * gsources)
                }
            }
            e[includes,,drop=nargs()==3L]   # drop = TRUE for loglik
        }
    }


    ### Calculates the two components of the integrated intensity function
    ### over [0;uppert] x W x K

    heIntTWK <- function (beta0, beta, gammapred, siafpars, tiafpars,
                          uppert = NULL) {}
    body(heIntTWK) <- as.call(c(as.name("{"),
        if (hash) { # endemic component
            expression(
                hIntTW <- .hIntTW(beta, uppert = uppert),
                .beta0 <- rep(if (nbeta0==0L) 0 else beta0, length.out = nTypes),
                fact <- sum(exp(.beta0)),
                hInt <- fact * hIntTW
            )
        } else { expression(hInt <- 0) },
        if (hase) { # epidemic component
            expression(
                siafInt <- do.call("..siafInt", .siafInt.args), # N-vector
                if (!is.null(uppert)) { # && isScalar(uppert) && t0 <= uppert && uppert < T
                    gIntUpper <- pmin(uppert-eventTimes, eps.t)
                    subtimeidx <- eventTimes < uppert
                    tiafIntSub <- .tiafInt(tiafpars,
                                           from = gIntLower[subtimeidx],
                                           to   = gIntUpper[subtimeidx],
                                           type = eventTypes[subtimeidx])
                    eInt <- sum(qSum[subtimeidx] * gammapred[subtimeidx] *
                                siafInt[subtimeidx] * tiafIntSub)
                } else {
                    tiafInt <- .tiafInt(tiafpars)
                    eInt <- sum(qSum * gammapred * siafInt * tiafInt)
                }
            )
        } else expression(eInt <- 0),
        expression(c(hInt, eInt))
    ))


    ### Calculates the log-likelihood

    loglik <- function (theta)
    {
        # Extract parameters from theta
        beta0    <- theta[seq_len(nbeta0)]
        beta     <- theta[nbeta0+seq_len(p)]
        gamma    <- theta[nbeta0+p+seq_len(q)]
        siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

        # dN part of the log-likelihood
        hEvents <- if (hash) .hEvents(beta0, beta) else 0
        eEvents <- if (hase) {
                gammapred <- drop(exp(mme %*% gamma)) # N-vector
                .eEvents(gammapred, siafpars, tiafpars) # Nin-vector! (only 'includes' here)
            } else 0
        lambdaEvents <- hEvents + eEvents  # Nin-vector
        llEvents <- sum(log(lambdaEvents))
        # here one might have got -Inf values in case of 0-intensity at an event time

        # lambda integral of the log-likelihood
        heInt <- heIntTWK(beta0, beta, gammapred, siafpars, tiafpars)   # !hase => missing(gammapred), but lazy evaluation omits an error in this case because heIntTWK doesn't ask for gammapred
        llInt <- sum(heInt)

        # Return the log-likelihood
        ll <- llEvents - llInt
        ll
    }


    ### Calculates the score vector

    score <- function (theta)
    {
        # Extract parameters from theta
        beta0    <- theta[seq_len(nbeta0)]
        beta     <- theta[nbeta0+seq_len(p)]
        gamma    <- theta[nbeta0+p+seq_len(q)]
        siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

        if (hase) {
            gammapred <- drop(exp(mme %*% gamma))  # N-vector
            hEvents <- if (hash) .hEvents(beta0, beta) else 0
            eEvents <- .eEvents(gammapred, siafpars, tiafpars) # Nin-vector! (only 'includes' here)
            lambdaEvents <- hEvents + eEvents  # Nin-vector
            siafInt <- do.call("..siafInt", .siafInt.args) # N-vector
            tiafInt <- .tiafInt(tiafpars) # N-vector
        }

        # score vector for beta
        hScore <- if (hash)
        {
            score_beta0 <- if (nbeta0 == 1L) local({ # global intercept
                sEvents <- if (hase) {
                        hEvents / lambdaEvents
                    } else rep.int(1, Nin)
                sEventsSum <- sum(sEvents)
                sInt <- nTypes*exp(beta0) * .hIntTW(beta)
                sEventsSum - sInt
            }) else if (nbeta0 > 1L) local({ # type-specific intercepts
                ind <- sapply(1:nTypes, function (type) eventTypes == type)
                #<- logical N x nTypes matrix
                sEvents <- if (hase) {
                        ind * hEvents / lambdaEvents
                    } else ind
                sEventsSum <- colSums(sEvents)
                sInt <- exp(beta0) * .hIntTW(beta)
                sEventsSum - sInt
            }) else numeric(0L) # i.e. nbeta0 == 0L

            score_beta <- if (p > 0L) local({
                sEvents <- if (hase) {
                        mmhEvents * hEvents / lambdaEvents
                    } else mmhEvents
                sEventsSum <- colSums(sEvents)
                fact <- if (nbeta0 > 1L) sum(exp(beta0)) else if (nbeta0 == 1L) nTypes*exp(beta0) else nTypes
                sInt <- fact * .hIntTW(beta, mmhGrid)
                sEventsSum - sInt
            }) else numeric(0L)

            c(score_beta0, score_beta)
        } else numeric(0L)

        # score vector for gamma, siafpars and tiafpars
        eScore <- if (hase)
        {
            score_gamma <- local({
                nom <- .eEvents(gammapred, siafpars, tiafpars, ncolsRes=q, score = mme)  # Ninxq matrix
                sEventsSum <- colSums(nom / lambdaEvents)
                sInt <- colSums(mme * (qSum * gammapred * siafInt * tiafInt))
                sEventsSum - sInt
            })

            score_siafpars <- if (hassiafpars && !fixedsiafpars) local({
                nom <- .eEvents(gammapred, siafpars, tiafpars,
                                ncolsRes=nsiafpars, f=siaf$deriv) # Nin x nsiafpars matrix
                sEventsSum <- colSums(nom / lambdaEvents)
                derivInt <- do.call(".siafDeriv", .siafDeriv.args) # N x nsiafpars matrix
                sInt <- colSums(derivInt * (qSum * gammapred * tiafInt))
                sEventsSum - sInt
            }) else numeric(nsiafpars) # if 'fixedsiafpars', this part is unused

            score_tiafpars <- if (hastiafpars && !fixedtiafpars) local({
                nom <- .eEvents(gammapred, siafpars, tiafpars,
                                ncolsRes=ntiafpars, g=tiaf$deriv) # Nin x ntiafpars matrix
                sEventsSum <- colSums(nom / lambdaEvents)
                derivIntUpper <- tiaf$Deriv(gIntUpper, tiafpars, eventTypes)
                derivIntLower <- tiaf$Deriv(gIntLower, tiafpars, eventTypes)
                derivInt <- derivIntUpper - derivIntLower
                sInt <- colSums(derivInt * (qSum * gammapred * siafInt))
                sEventsSum - sInt
            }) else numeric(ntiafpars) # if 'fixedtiafpars', this part is unused

            c(score_gamma, score_siafpars, score_tiafpars)
        } else numeric(0L)

        # return the score vector
        scorevec <- c(hScore, eScore)
        scorevec
    }


    ### Estimates the expected Fisher information matrix
    ### by the "optional variation process" (Martinussen & Scheike, p. 64),
    ### or see Rathbun (1996, equation (4.7))

    fisherinfo <- function (theta)
    {
        # Extract parameters from theta
        beta0    <- theta[seq_len(nbeta0)]
        beta     <- theta[nbeta0+seq_len(p)]
        gamma    <- theta[nbeta0+p+seq_len(q)]
        siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

        # only events (intdN) part of the score function needed
        zeromatrix <- matrix(0, Nin, 0)

        if (hase) {
            gammapred <- drop(exp(mme %*% gamma))  # N-vector
            hEvents <- if (hash) .hEvents(beta0, beta) else 0
            eEvents <- .eEvents(gammapred, siafpars, tiafpars) # Nin-vector! (only 'includes' here)
            lambdaEvents <- hEvents + eEvents  # Nin-vector
        }

        # for beta
        hScoreEvents <- if (hash) {
            scoreEvents_beta0 <- if (nbeta0 > 1L) local({ # type-specific intercepts
                ind <- sapply(1:nTypes, function (type) eventTypes == type) # logical NxnTypes matrix
                if (hase) {
                    ind * hEvents / lambdaEvents
                } else ind
            }) else if (nbeta0 == 1L) { # global intercept
                if (hase) {
                    hEvents / lambdaEvents
                } else matrix(1, Nin, 1L)
            } else zeromatrix

            scoreEvents_beta <- if (p > 0L) {
                if (hase) {
                    mmhEvents * hEvents / lambdaEvents
                } else mmhEvents   # Ninxp matrix
            } else zeromatrix

            cbind(scoreEvents_beta0, scoreEvents_beta)
        } else zeromatrix

        # for gamma, siafpars and tiafpars
        eScoreEvents <- if (hase)
        {
            scoreEvents_gamma_nom <-
                .eEvents(gammapred, siafpars, tiafpars, ncolsRes = q, score = mme)  # Ninxq matrix

            scoreEvents_siafpars_nom <- if (hassiafpars) {
                .eEvents(gammapred, siafpars, tiafpars, ncolsRes = nsiafpars, f = siaf$deriv)  # Ninxnsiafpars matrix
            } else zeromatrix

            scoreEvents_tiafpars_nom <- if (hastiafpars) {
                .eEvents(gammapred, siafpars, tiafpars, ncolsRes = ntiafpars, g = tiaf$deriv)  # Ninxntiafpars matrix
            } else zeromatrix

            eScoreEvents_nom <- cbind(scoreEvents_gamma_nom, scoreEvents_siafpars_nom, scoreEvents_tiafpars_nom)
            eScoreEvents_nom / lambdaEvents
        } else zeromatrix

        scoreEvents <- cbind(hScoreEvents, eScoreEvents)

        # Build the optional variation process (Martinussen & Scheike, p64)
        info <- matrix(0, nrow = npars, ncol = npars,
                    dimnames = list(names(theta), names(theta)))
        for (i in 1:Nin) {
            x <- scoreEvents[i,,drop=FALSE]  # single-ROW matrix
            info <- info + crossprod(x) # t(x) %*% x
        }

        # Return the estimated Fisher information matrix
        info
    }


    ### Calculates the partial log-likelihood for continuous space
    ### (Diggle et al., 2009)

    partialloglik <- function (theta)
    {
        # Extract parameters from theta
        beta0    <- theta[seq_len(nbeta0)]
        beta     <- theta[nbeta0+seq_len(p)]
        gamma    <- theta[nbeta0+p+seq_len(q)]
        siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

        # calculcate the observed intensities
        hEvents <- if (hash) .hEvents(beta0, beta) else 0
        eEvents <- if (hase) {
                gammapred <- drop(exp(mme %*% gamma))  # N-vector
                .eEvents(gammapred, siafpars, tiafpars)  # Nin-vector! (only 'includes' here)
            } else 0
        lambdaEvents <- hEvents + eEvents  # Nin-vector

        # calculate integral of lambda(t_i, s, kappa) over at-risk set = (observation region x types)
        hInts <- if (hash) { # endemic component
                etahGrid <- if (p > 0L) drop(mmhGrid %*% beta) else numeric(nrow(mmhGrid))
                if (!is.null(offsetGrid)) etahGrid <- offsetGrid + etahGrid
                hGrid <- exp(etahGrid)
                # integral over W and types for each time block in mfhGrid
                fact <- if (nbeta0 > 1L) sum(exp(beta0)) else if (nbeta0 == 1L) nTypes*exp(beta0) else nTypes
                hInt_blocks <- fact * tapply(hGrid*ds, gridBlocks, sum, simplify = TRUE)
                .idx <- match(eventBlocks[includes], names(hInt_blocks))
                unname(hInt_blocks[.idx])   # Nin-vector
            } else 0
        eInts <- if (hase) { # epidemic component
                siafInt <- do.call(".siafInt", .siafInt.args) # N-vector
                gs <- gammapred * siafInt # N-vector
                sapply(includes, function (i) {
                    timeSources <- determineSources(i, eventTimes, removalTimes,
                        0, Inf, NULL)
                    nSources <- length(timeSources)
                    if (nSources == 0L) 0 else {
                        repi <- rep.int(i, nSources)
                        tdiff <- eventTimes[repi] - eventTimes[timeSources]
                        gsources <- tiaf$g(tdiff, tiafpars, eventTypes[timeSources])
                        sum(qSum[timeSources] * gs[timeSources] * gsources)
                    }
                })   # Nin-vector
            } else 0
        lambdaEventsIntW <- hInts + eInts   # Nin-vector

        # Calculate and return the partial log-likelihood
        p <- lambdaEvents / lambdaEventsIntW   # Nin-vector
        pll <- sum(log(p))
        pll
    }






    ################################
    ### Prepare for optimization ###
    ################################


    ll <- if (partial) partialloglik else loglik
    functions <- list(ll = ll,
                      sc = if (useScore) score else NULL,
                      fi = if (useScore) fisherinfo else NULL)

    
    ### Include check for validity of siafpars and tiafpars ('validpars') in ll

    if (!is.null(siaf$validpars)) {
        body(ll) <- as.call(append(as.list(body(ll)),
            as.list(expression(
                if (hassiafpars && !siaf$validpars(siafpars)) {
                    if (optimArgs$control$trace > 0L)
                        cat("(invalid 'siafpars' in loglik)\n")
                    return(-Inf)
                }
                )),
            after = grep("^siafpars <-", body(ll))))
    }

    if (!is.null(tiaf$validpars)) {
        body(ll) <- as.call(append(as.list(body(ll)),
            as.list(expression(
                if (hastiafpars && !tiaf$validpars(tiafpars)) {
                    if (optimArgs$control$trace > 0L)
                        cat("(invalid 'tiafpars' in loglik)\n")
                    return(-Inf)
                }
                )),
            after = grep("^tiafpars <-", body(ll))))
    }


    ### Check that optim.args is a list or NULL

    if (missing(optim.args) || (!is.list(optim.args) && !is.null(optim.args))) {
        stop("'optim.args' must be a list or NULL")
    }

    if (is.null(optim.args)) {          # no optimisation requested
        setting <- functions
        on.exit(rm(setting), add = TRUE)
        # Append model information
        setting$npars <- c(nbeta0 = nbeta0, p = p,
                           q = q, nsiafpars = nsiafpars, ntiafpars = ntiafpars)
        setting$qmatrix <- qmatrix   # -> information about nTypes and typeNames
        setting$formula <- list(endemic = formula(endemic), epidemic = formula(epidemic),
                            siaf = siaf, tiaf = tiaf)
        setting$call <- cl
        # Return settings
        message("optimization skipped (returning functions in data environment)")
        return(setting)
    }

    
    ### Check start value for theta

    if (is.null(optim.args[["par"]])) {
        stop("start values of parameters must be specified as 'optim.args$par'")
    }
    if (!is.vector(optim.args$par, mode="numeric")) {
        stop("'optim.args$par' must be a numeric vector")
    }
    if (length(optim.args$par) != npars) {
        stop(gettextf(paste("'optim.args$par' (%d) does not have the same",
             "length as the number of unknown parameters (%d)"),
             length(optim.args$par), npars))
    }
    
    # Set names for theta
    names(optim.args$par) <- c(
        if (nbeta0 > 1L) {
            paste0("h.type",typeNames)
        } else if (nbeta0 == 1L) "h.(Intercept)",
        if (p > 0L) paste("h", colnames(mmhEvents), sep = "."),
        if (hase) paste("e", colnames(mme), sep = "."),
        if (hassiafpars) paste("e.siaf",1:nsiafpars,sep="."),
        if (hastiafpars) paste("e.tiaf",1:ntiafpars,sep=".")
    )


    ### Fixed parameters during optimization

    fixed <- optim.args[["fixed"]]
    optim.args[["fixed"]] <- NULL
    if (is.null(fixed)) fixed <- integer(0L) else stopifnot(is.vector(fixed))
    whichfixed <- if (is.numeric(fixed)) {
        stopifnot(fixed %in% 1:npars)
        fixed
    } else if (is.character(fixed)) {
        if (any(!fixed %in% names(optim.args$par)))
            stop("'optim.args$fixed' must be a subset of: ",
                 paste0("\"",names(optim.args$par), "\"", collapse=", "))
        fixed
    } else if (is.logical(fixed)) {
        stopifnot(length(fixed) == npars)
        which(fixed)
    } else stop("'optim.args$fixed' must be a numeric, character or logical vector")
    fixed <- logical(npars)   # FALSE
    names(fixed) <- names(optim.args$par)
    fixed[whichfixed] <- TRUE
    fixedsiafpars <- hassiafpars && all(fixed[paste("e.siaf", 1:nsiafpars, sep=".")])
    fixedtiafpars <- hastiafpars && all(fixed[paste("e.tiaf", 1:ntiafpars, sep=".")])

    ## in the end, we set fixed[st]iafpars to FALSE (for free posteriori evaluations)
    on.exit(fixedsiafpars <- fixedtiafpars <- FALSE, add = TRUE)


    ### Define negative log-likelihood (score, hessian) for minimization
    ### as a function of the non-fixed parameters
    
    negll <- ll
    body(negll)[[length(body(negll))]] <-
        call("-", body(negll)[[length(body(negll))]])
    negsc <- if (useScore) {
        negsc <- score
        body(negsc)[[length(body(negsc))]] <-
            call("-", body(negsc)[[length(body(negsc))]])
        negsc
    } else NULL
    neghess <- if (useScore) fisherinfo else NULL

    if (any(fixed)) {
        ## modify negll, negsc and neghess for subvector optimization
        initpars <- optim.args$par
        optim.args$par <- initpars[!fixed]
        if (all(fixed)) {
            cat("\nno numerical likelihood optimization, all parameters fixed:\n")
        } else cat("\nfixed parameters during optimization:\n")
        print(initpars[fixed])
        tmpexpr <- expression(
            initpars[!fixed] <- theta,
            theta <- initpars
            )
        body(negll) <- as.call(append(as.list(body(negll)), as.list(tmpexpr), 1))
        if (useScore) {
            body(negsc) <- as.call(append(as.list(body(negsc)), as.list(tmpexpr), 1))
            body(neghess) <- as.call(append(as.list(body(neghess)), as.list(tmpexpr), 1))
            # return non-fixed sub-vector / sub-matrix only
            body(negsc)[[length(body(negsc))]] <-
                call("[", body(negsc)[[length(body(negsc))]], quote(!fixed))
            body(neghess)[[length(body(neghess))]] <-
                call("[", body(neghess)[[length(body(neghess))]],
                     quote(!fixed), quote(!fixed), drop=FALSE)
        }

        ## if siafpars or tiafpars are fixed, pre-evaluate integrals    
        if (fixedsiafpars) {
            cat("pre-evaluating 'siaf' integrals with fixed parameters ...\n")
            .siafInt.args[[1]] <- initpars[paste("e.siaf", 1:nsiafpars, sep=".")]
            siafInt <- do.call(".siafInt", .siafInt.args)
            ## re-define .siafInt such that it just returns the pre-evaluated integrals
            .siafInt.orig <- .siafInt
            body(.siafInt) <- expression(siafInt)
            ## restore the original function at the end
            on.exit({
                .siafInt <- .siafInt.orig
                rm(.siafInt.orig)
            }, add=TRUE)
        }
        if (fixedtiafpars) {
            cat("pre-evaluating 'tiaf' integrals with fixed parameters ...\n")
            tiafInt <- .tiafInt(initpars[paste("e.tiaf", 1:ntiafpars, sep=".")])
            ## re-define .tiafInt such that it just returns the pre-evaluated
            ## integrals if called with the default arguments
            .tiafInt.orig <- .tiafInt
            body(.tiafInt) <- expression(
                if (nargs() == 1L) tiafInt else
                .tiafInt.orig(tiafpars, from, to, type, G)
            )
            ## restore the original function at the end
            on.exit({
                .tiafInt <- .tiafInt.orig
                rm(.tiafInt.orig)
            }, add=TRUE)
        }
    }






    if (any(!fixed)) {   ####################
                         ### Optimization ###
                         ####################

        ## Configure the optim procedure (check optim.args)

        # default arguments
        optimArgs <- alist(par =, fn = negll, gr = negsc,
                           method = if (partial) "Nelder-Mead" else "nlminb",
                           lower = -Inf, upper = Inf,
                           control = list(), hessian = TRUE)
        # user arguments
        namesOptimArgs <- names(optimArgs)
        namesOptimUser <- names(optim.args)
        optimValid <- namesOptimUser %in% namesOptimArgs
        optimArgs[namesOptimUser[optimValid]] <- optim.args[optimValid]
        if (any(!optimValid)) {
            warning("unknown names in optim.args: ",
                    paste(namesOptimUser[!optimValid], collapse = ", "))
        }
        doHessian <- eval(optimArgs$hessian)
        optimMethod <- eval(optimArgs$method)


        ## Call 'optim', 'nlminb', or 'nlm' with the above arguments

        cat("\nminimizing the negative", if (partial) "partial", "log-likelihood",
            "using", if (optimMethod %in% c("nlm", "nlminb"))
            paste0("'",optimMethod,"()'") else {
                paste0("'optim()'s \"", optimMethod, "\"")
            }, "...\n")
        cat("initial parameters:\n")
        print(optimArgs$par)
        optimRes1 <- if (optimMethod == "nlminb") {
            nlminbControl <- control2nlminb(optimArgs$control,
                                            defaults = list(trace=5L, rel.tol=1e-6))
            ## sqrt(.Machine$double.eps) is the default reltol used in optim,
            ## which usually equals about 1.49e-08.
            ## The default rel.tol of nlminb (1e-10) seems too small
            ## (nlminb often does not finish despite no "relevant" change in loglik).
            ## I therefore use 1e-6, which is also the default in package nlme
            ## (see 'lmeControl').
            if (nlminbControl$trace > 0L) {
                cat("negative log-likelihood and parameters ")
                if (nlminbControl$trace == 1L) cat("in each iteration") else {
                    cat("every", nlminbControl$trace, "iterations") }
                cat(":\n")
            }
            nlminbRes <- nlminb(start = optimArgs$par, objective = negll,
                                gradient = negsc,
                                hessian = if (doHessian) neghess else NULL,
                                control = nlminbControl,
                                lower = optimArgs$lower, upper = optimArgs$upper)
            nlminbRes$value <- -nlminbRes$objective
            nlminbRes$counts <- nlminbRes$evaluations
            nlminbRes
        } else if (optimMethod == "nlm") {
            nlmObjective <- function (theta) {
                value <- negll(theta)
                grad <- negsc(theta)
                #hess <- neghess(theta)
                structure(value, gradient = grad)#, hessian = hess)
            }
            nlmControl <- optimArgs$control
            if (is.null(nlmControl[["print.level"]])) {
                nlmControl$print.level <- min(nlmControl$trace, 2L)
            }
            nlmControl$trace <- nlmControl$REPORT <- NULL
            if (is.null(nlmControl[["iterlim"]])) {
                nlmControl$iterlim <- nlmControl$maxit
            }
            nlmControl$maxit <- NULL
            nlmControl$check.analyticals <- FALSE
            ##<- we use the negative _expected_ Fisher information as the Hessian,
            ##   which is of course different from the true Hessian (=neg. obs. Fisher info)
            nlmRes <- do.call("nlm", c(alist(f = nlmObjective, p = optimArgs$par,
                                             hessian = doHessian),
                                             nlmControl))
            names(nlmRes)[names(nlmRes) == "estimate"] <- "par"
            nlmRes$value <- -nlmRes$minimum
            nlmRes$counts <- rep.int(nlmRes$iterations, 2L)
            nlmRes$convergence <- if (nlmRes$code %in% 1:2) 0L else nlmRes$code
            nlmRes
        } else { # use optim()
            optimArgs$control <- modifyList(list(trace=1L, REPORT=5L),
                                            optimArgs$control)
            if (finetune) optimArgs$hessian <- FALSE
            res <- do.call("optim", optimArgs)
            res$value <- -res$value
            res
        }


        ## Optional fine-tuning of ML estimates by robust Nelder-Mead

        optimRes <- if (finetune) {
            cat("\nMLE from first optimization:\n")
            print(optimRes1$par)
            cat("loglik(MLE) =", optimRes1$value, "\n")
            cat("\nfine-tuning MLE using Nelder-Mead optimization...\n")
            optimArgs$par <- optimRes1$par
            optimArgs$method <- "Nelder-Mead"
            optimArgs$hessian <- doHessian
            optimArgs$control <- modifyList(list(trace=1L), optimArgs$control)
            nmRes <- do.call("optim", optimArgs)
            nmRes$value <- -nmRes$value
            nmRes$counts[2L] <- 0L   # 0 gradient evaluations (replace NA for addition below)
            nmRes
        } else optimRes1

        if (optimRes$convergence != 0) {
            cat("\nWARNING: OPTIMIZATION ROUTINE DID NOT CONVERGE",
                if (finetune || optimMethod != "nlminb")
                paste0("(code ", optimRes$convergence, ")"),
                "!\n")
            if (!is.null(optimRes$message) && nzchar(optimRes$message)) {
                cat("MESSAGE: \"", optimRes$message, "\"\n", sep="")
            }
        }

        cat("\n", if (finetune) "final ", "MLE:\n", sep = "")
        print(optimRes$par)
        cat("loglik(MLE) =", optimRes$value, "\n")

    }






    ##############
    ### Return ###
    ##############


    ## ### Issue collected polyCub.midpoint warnings
    
    ## sapply(unique(.POLYCUB.WARNINGS), warning)


    ### Set up list object to be returned

    fit <- list(
           coefficients = if (any(fixed)) {
               if (all(fixed)) initpars else
               unlist(modifyList(as.list(initpars), as.list(optimRes$par)))
           } else optimRes$par,
           loglik = structure(if (all(fixed)) ll(initpars) else optimRes$value,
                              partial = partial),
           counts = if (all(fixed)) c("function"=1L, "gradient"=0L) else {
               optimRes1$counts + if (finetune) optimRes$counts else c(0L, 0L)
           },
           converged = all(fixed) || (optimRes$convergence == 0)
           )


    ### Add Fisher information matrices

    # estimation of the expected Fisher information matrix
    if (useScore) fit$fisherinfo <- fisherinfo(fit$coefficients)

    # If requested, add observed fisher info (= negative hessian at maximum)
    if (any(!fixed) && !is.null(optimRes$hessian)) {
        fit$fisherinfo.observed <- optimRes$hessian
        ## no "-" here because we optimized the negative log-likelihood
    }


    ### Add fitted intensity values and integrated intensities at events

    # final coefficients
    theta    <- fit$coefficients
    beta0    <- theta[seq_len(nbeta0)]
    beta     <- theta[nbeta0+seq_len(p)]
    gamma    <- theta[nbeta0+p+seq_len(q)]
    siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
    tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

    # final siaf and tiaf integrals over influence regions / periods
    # and final gammapred (also used by intensity.twinstim)
    if (hase) {
        gammapred <- drop(exp(mme %*% gamma)) # N-vector
        if (!fixedsiafpars) siafInt <- do.call("..siafInt", .siafInt.args)
        if (!fixedtiafpars) tiafInt <- .tiafInt(tiafpars)
    }
    
    # fitted intensities
    hEvents <- if (hash) .hEvents(unname(beta0), beta) else rep.int(0, Nin)
    eEvents <- if (hase) {
            .eEvents(gammapred, siafpars, tiafpars) # Nin-vector! (only 'includes' here)
        } else rep.int(0, Nin)
    fit$fitted <- hEvents + eEvents   # = lambdaEvents  # Nin-vector
    fit$fittedComponents <- cbind(h = hEvents, e = eEvents)
    rm(hEvents, eEvents)
    
    # calculate cumulative ground intensities at event times
    # Note: this function is also used by residuals.twinstim
    LambdagEvents <- function (cumCIF.pb = TRUE)
    {
        if (hase) {
            ## tiny hack such that siaf integrals are not evaluated N-fold
            .siafInt.orig <- .siafInt
            body(.siafInt) <<- expression(siafInt)
            on.exit(.siafInt <<- .siafInt.orig)
        }
        heIntEvents <- matrix(NA_real_, Nin, 2L)
        if (cumCIF.pb) pb <- txtProgressBar(min=0, max=Nin, initial=0, style=3)
        for (i in 1:Nin) {
            heIntEvents[i,] <- heIntTWK(beta0, beta, gammapred, siafpars,
                                        tiafpars, uppert = eventTimes[includes[i]])
            if (cumCIF.pb) setTxtProgressBar(pb, i)
        }
        if (cumCIF.pb) close(pb)
        LambdagEvents <- rowSums(heIntEvents)
        names(LambdagEvents) <- rownames(mmhEvents) # this is NULL if only h.intercept
                                                    # (cf. bug report #14992)
        if (is.null(names(LambdagEvents))) names(LambdagEvents) <- rownames(mme)
        LambdagEvents
    }
    if (cumCIF) {
        cat("\nCalculating the fitted cumulative intensities at events...\n")
        fit$tau <- LambdagEvents(cumCIF.pb)
    }

    # calculate observed R0's: mu_j = spatio-temporal integral of e_j(t,s) over
    # the observation domain (t0;T] x W (not whole R+ x R^2)
    fit$R0 <- if (hase) qSum * gammapred * siafInt * tiafInt else rep.int(0, N)
    names(fit$R0) <- rownames(mfe)

    
    ### Append model information

    fit$npars <- c(nbeta0 = nbeta0, p = p,
                   q = q, nsiafpars = nsiafpars, ntiafpars = ntiafpars)
    fit$qmatrix <- qmatrix   # -> information about nTypes and typeNames
    fit$bbox <- bbox(data$W)            # for completeness and for iafplot
    fit$timeRange <- c(t0, T)           # for simulate.twinstim's defaults
    if (!model) {
        # Link formulae to the global environment such that the evaluation
        # environment will be dropped at the end
        environment(epidemic) <- environment(endemic) <- .GlobalEnv
    }
    # if typeSpecificEndemicIntercept, re-add this to the endemic formula
    fit$formula <- list(
                   endemic = if (typeSpecificEndemicIntercept) {
                       update.formula(formula(endemic), ~ (1|type) + .)
                   } else formula(endemic),
                   epidemic = formula(epidemic),
                   siaf = siaf, tiaf = tiaf
                   )
    fit$control.siaf <- control.siaf

    
    ### Append optimizer configuration
    
    fit$optim.args <- optim.args
    if (model) {
        fit$functions <- functions
        environment(fit) <- environment()
    }


    ### Return object of class "twinstim"

    cat("\nDone.\n")
    fit$call <- cl
    fit$runtime <- proc.time()[[3]] - ptm
    class(fit) <- "twinstim"
    return(fit)

}
