\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
%% some pre-defined commands: \R, \code, \acronym, \url, \file, \pkg

\name{NEWS}
\title{News for Package 'surveillance'}
\encoding{latin1}


\section{Changes in surveillance version 1.10-0 (2015-11-04)}{

\subsection{NEW FEATURES}{
  \itemize{
    \item Calibration tests for count data (Wei and Held, 2014, Test)
    are now implemented and available as \code{calibrationTest()}.
    In addition to a default method taking pure counts and predictive
    means and dispersion parameters, there are convenient methods for
    \code{"hhh4"} and \code{"oneStepAhead"} objects.

    \item Shared overdispersion across units in negative binomial
    \code{hhh4()} time series models (by specifying a factor variable
    as the \code{family} argument).

    \item \code{scores()} and \code{pit()} are now generic and have convenient
    methods for \code{"oneStepAhead"} predictions and \code{"hhh4"} fits.

    \item The initial values used for model updates during the
    \code{oneStepAhead()} procedure can now be specified directly
    through the \code{which.start} argument (as an alternative to the
    previous options \code{"current"} and \code{"final"}).
    
    \item \code{plotHHH4_fitted()} (and \code{plotHHH4_fitted1()})
    gained an option \code{decompose} to plot the contributions from
    each single unit (and the endemic part) instead of the default
    endemic + AR + neighbours decomposition.
    Furthermore, a formatted time axis similar to \code{stsplot_time1()}
    can now be enabled via the new argument \code{xaxis}.

    \item The new \code{plot} \code{type} \code{"maps"} for \code{"hhh4"}
    fits shows maps of the fitted mean components averaged over time.

    \item New \code{plot}-method for simulations from \code{"hhh4"}
    models (using \code{simulate.hhh4(..., simplify = TRUE)}, which now
    has a dedicated class: \code{"hhh4sims"}) to show the final
    size distribution or the simulated time series
    (possibly stratified by groups of units).
    There is also a new \code{scores}-method to compute proper scoring
    rules based on such simulations.

    \item The argument \code{idx2Exp} of \code{coef.hhh4()} may now be
    conveniently set to \code{TRUE} to exp-transform all coefficients.

    \item Added a \code{coeflist()}-method for \code{"hhh4"} fits.

    \item The generator function \code{sts()} can now be used to
    initialize objects of class \code{"sts"} (instead of writing
    \code{new("sts", ...)}).

    \item Additional arguments of \code{layout.scalebar()} now allow to
    change the style of the labels.

    \item A pre-computed distance matrix \code{D} can now be used as
    input for the \code{as.epidata()} converter -- offering an alternative
    to the default Euclidean distance based on the individuals coordinates.
    (Request of George Wood to support \code{twinSIR} models on networks.)
  }
}

\subsection{MINOR CHANGES}{
  \itemize{
    \item The first argument of \code{scores()} is now called \code{x}
    instead of \code{object} (for consistency with \code{calibrationTest()}).

    \item The result of \code{oneStepAhead()} now has the dedicated
    class attribute \code{"oneStepAhead"} (previously was just a list).

    \item Changed interpretation of the \code{col} argument of
    \code{plotHHH4_fitted()} and \code{plotHHH4_fitted1()} (moved color
    of \dQuote{observed} to separate argument \code{pt.col} and reversed
    remaining colors). The old \code{col} specification as a vector of
    length 4 still works (catched internally) but is undocumented.

    \item The \code{epoch} slot of class \code{"sts"} is now initialized to
    \code{1:nrow(observed)} by default and thus no longer needs to be
    explicitly set when creating a \code{new("sts", ...)} for this
    standard case.

    \item Initialization of \code{new("sts", ...)} now supports the
    argument \code{frequency} (for consistency with \code{ts()}).
    Note that \code{freq} still works (via partial argument matching)
    and that the corresponding \code{"sts"} slot is still called \code{freq}.

    \item If \code{missing(legend.opts)} in \code{stsplot_time1()}, the
    default legend will only be produced if the \code{"sts"} object
    contains information on outbreaks, alarms, or upperbounds.

    \item The default \code{summary()} of a \code{"twinstim"} fit is more
    concise since it no longer includes the number of log-likelihood and
    score function evaluations and the elapsed time during model fitting.
    Set the new \code{runtime} argument of \code{summary.twinstim()} to
    \code{TRUE} to add this information to the summary as before.

    \item The \code{animate}-method for \code{"sts"} objects gained an
    argument \code{draw} (to disable the default instantaneous plotting)
    and now invisibly returns the sequential plot objects (of class
    \code{"gtable"} or \code{"trellis"}) in a list for post-processing.

    \item The flexible time axis configurations for \code{"sts"} plots
    introduced in version 1.8-0 now also work for classical \code{"sts"}
    objects with integer epochs and standard frequencies
    (try \code{plot(..., epochsAsDate = TRUE)}).

    \item \code{stsplot_time()} initiates \code{par} settings only
    if the \code{par.list} argument is a list.

    \item The new \code{all.equal()} method for class \code{"hhh4"}
    compares two fits ignoring their \code{"runtime"} and \code{"call"}
    elements (at least).
  }
}

\subsection{BUG FIXES}{
  \itemize{
    \item Fixed a bug in \code{algo.bayes}, where an alarm was already
    sounded if the current observation was equal to the quantile of the
    predictive posterior. This was changed in order to get \eqn{alarm_t
    = I(obs_t > quantile_t)} which is consistent with the use in
    \code{boda} and \code{bodaDelay}.

    \item Fixed bug in \code{algo.outbreakP} causing a halt in the
    computations of \code{value="cases"} when
    \code{calc.outbreakP.statistic} returned \code{NaN}. Now, a
    \code{NaN} is returned.

    \item \code{wrap.algo} argument \code{control.hook} used
    \code{control} argument defined outside it's scope (and not the one
    provided to the function). It is now added as additional 2nd argument
    to the \code{control.hook} function.

    \item \code{stsplot_time()} did not account for the optional
    \code{units} argument for multivariate \code{"sts"} objects
    when choosing a suitable value for \code{par("mfrow")}.

    \item \code{hhh4()} could have used a function \code{dpois()} or
    \code{dnbinom()} from the global environment instead of the
    respective function from package \pkg{stats}.

    \item The default time variable \code{t} created as part of the
    \code{data} argument in \code{hhh4()} was incompatible with
    \code{"sts"} objects having \code{epochAsDate=TRUE}.

    \item A consistency check in \code{as.epidata.default()} failed for
    SI-type data (and, more generally, for all data which ended with an
    I-event in the last time block). Spotted by George Wood.
  }
}

}



\section{Changes in surveillance version 1.9-1 (2015-06-12)}{

  \itemize{
    \item This is a quick patch release to make the test suite run
    smoothly on CRAN's Windows and Solaris Sparc systems.

    \item The new \code{hhh4()} option to scale neighbourhood weights
    did not work for parametric weights with more than one parameter
    if \code{normalize=FALSE}.
  }

}



\section{Changes in surveillance version 1.9-0 (2015-06-09)}{

\subsection{NEW FEATURES}{
  \itemize{
    \item New functions and data for Bayesian outbreak detection in the
    presence of reporting delays (Salmon et al., 2015):
    \code{bodaDelay()}, \code{sts_observation()}, and \code{sts_creation()}.

    \item New functions implementing tests for space-time interaction:
    \itemize{
      \item \code{knox()} supports both the Poisson approximation and a
      Monte Carlo permutation approach to determine the p-value,
      \item \code{stKtest()} wraps space-time K-function methods from
      package \CRANpkg{splancs} for use with \code{"epidataCS"},
      \item and \code{epitest()} for \code{twinstim} models
      (makes use of the new auxiliary function \code{simpleR0()}).
    }

    \item New function \code{plapply()}: a parallel and verbose version
    of \code{lapply()} wrapping around both \code{mclapply()} and
    \code{parLapply()} of package \pkg{parallel}.

    \item New converter \code{as.xts.sts()} to transform \code{"sts"}
    objects to the quasi standard \code{"xts"} class, e.g., to make use
    of package \CRANpkg{dygraphs} for interactive time series plots.

    \item New options for scaling and normalization of neighbourhood
    weights in \code{hhh4()} models.

    \item New auxiliary function \code{layout.scalebar()} for use as part
    of \code{sp.layout} in \code{spplot()} or in the traditional
    graphics system.
  }

  \subsection{New features for \code{"epidataCS"}}{
  \itemize{
    \item New argument \code{by} for \code{plot.epidataCS()}, which
    defines a stratifying variable for the events (default is the event
    type as before). It can also be set to \code{NULL} to make the plot
    not distinguish between event types.

    \item The spatial plot of \code{"epidataCS"} gained the arguments
    \code{tiles}, \code{pop} and \code{sp.layout}, and can now produce
    an \code{spplot()} with the tile-specific population levels behind
    the point pattern.

    \item New function \code{permute.epidataCS()} to randomly permute
    time points or locations of the events (holding other marks fixed).
  }
  }

  \subsection{New features for \code{twinstim()}}{
  \itemize{
    \item New S3-generic \code{coeflist()} to list model coefficients by
    component. It currently has a default method and one for
    \code{"twinstim"} and \code{"simEpidataCS"}.

    \item New argument \code{newcoef} for \code{simulate.twinstim()} to
    customize the model parameters used for the simulation.

    \item New argument \code{epilink} for \code{twinstim()}, offering
    experimental support for an identity link for the epidemic
    predictor. The default remains \code{epilink = "log"}.

    \item Simulation from \code{"twinstim"} models and generation of
    \code{"epidataCS"} is slightly faster now (faster \pkg{spatstat}
    functions are used to determine the distance of events to the border).

    \item New option \code{scaled = "standardized"} in \code{iafplot()}
    to plot \eqn{f(x) / f(0)} or \eqn{g(t) / g(0)}, respectively.
  }
  }
}

\subsection{MINOR CHANGES}{
  \itemize{
    \item Initial data processing in \code{twinstim()} is faster
    since event sources are only re-determined if there is effective
    need for an update (due to subsetting or a change of
    \code{qmatrix}).

    \item \code{formatPval()} disables \code{scientific} notation by default.

    \item The \code{"time"} plot for \code{"epidataCS"} uses the
    temporal grid points as the default histogram \code{breaks}.

    \item The special \code{fe()} function which sets up fixed effects
    in \code{hhh4()} models gained an argument \code{unitSpecific} as a
    convenient shortcut for \code{which = rep(TRUE, nUnits)}.

    \item The convenient \code{plot} option of \code{permutationTest()}
    uses \CRANpkg{MASS}::\code{truehist()} instead of \code{hist()} and
    accepts graphical parameters to customize the histogram.
  }
}

\subsection{BUG FIXES}{
  \itemize{
    \item The \code{bodaFit} function did not draw samples from the
    joint posterior. Instead draws were from the respective posterior
    marginals. A new argument \code{samplingMethod} is now introduced
    defaulting to the proper 'joint'. For backwards compatibility use
    the value 'marginal'.

    \item The functions \code{as.epidataCS()} and \code{simEpidataCS()}
    could throw inappropriate warnings when checking polygon areas
    (only if \code{W} or \code{tiles}, respectively, contained holes).

    \item Non-convergent endemic-only \code{twinstim} models
    produced an error (spotted by Bing Zhang).

    \item The \code{"owin"}-method of \code{intersectPolyCircle} could
    have returned a rectangle-type \code{"owin"} instead of a polygon.

    \item An error occurred in \code{twinstim()} if \code{finetune=TRUE}
    or choosing \code{optim()} instead of the default \code{nlminb()}
    optimizer without supplying a \code{control} list in \code{optim.args}.

    \item The \code{"time"} plot for \code{"epidataCS"} did not
    necessarily use the same histogram \code{breaks} for all strata.

    \item Specifying a step function of interaction via a numeric vector
    of knots did not work in \code{twinstim()}.

    \item \code{plot.hhh4()} did not support an unnamed \code{type}
    argument such as \code{plot(x, "season")}.

    \item \code{simEpidataCS()} did not work if \code{t0} was in the
    last block of \code{stgrid} (thus it did not work for single-cell
    grids), and mislabeled the \code{start} column copied to
    \code{events} if there were no covariates in \code{stgrid}.

    \item Evaluating \code{intensity.twinstim()$hFUN()} at time points
    before \code{t0} was an error. The function now returns
    \code{NA_real_} as for time points beyond \code{T}.

    \item Truncated, normalized power-law weights for \code{hhh4()}
    models, i.e., \code{W_powerlaw(maxlag = M, normalize = TRUE)}
    with \code{M < max(neighbourhood(stsObj))}, had wrong derivatives
    and thus failed to converge.

    \item \code{update.hhh4(..., use.estimates = TRUE)} did not
    use the estimated weight function parameters as initial values for
    the new fit. It does so now iff the weight function
    \code{ne$weights} is left unchanged.
  }
}

}



\section{Changes in surveillance version 1.8-3 (2015-01-05)}{

  \itemize{
    \item Accommodate a new note given by R-devel checks, and set the new
    INLA additional repository in the \file{DESCRIPTION} file.

    \item Made \code{linelist2sts()} work for quarters by adding extra
    \code{"\%q"} formatting in \code{formatDate()}.
  }

}



\section{Changes in surveillance version 1.8-2 (2014-12-16)}{

\subsection{MINOR CHANGES related to \code{hhh4}}{
  \itemize{
    \item In the coefficient vector resulting from a \code{hhh4} fit,
    random intercepts are now named.

    \item Parameter \code{start} values in \code{hhh4()} are now
    matched by name but need not be complete in that case (default
    initial values are used for unspecified parameters).

    \item The \code{update.hhh4()}-method now by default does
    \code{use.estimates} from the previous fit. This reduces the
    number of iterations during model fitting but may lead to slightly
    different parameter estimates (within a tolerance of \code{1e-5}).
    Setting \code{use.estimates = FALSE} means to re-use the previous
    start specification.
  }
}

\subsection{MINOR CHANGES related to the \code{"sts"}-class}{
  \itemize{
    \item For univariate \code{"sts"} objects, the (meaningless)
    \dQuote{head of neighbourhood} is no longer \code{show}n.

    \item The \code{"sts"} class now has a \code{dimnames}-method
    instead of a \code{colnames}-method. Furthermore, the redundant
    \code{nrow} and \code{ncol} methods have been removed
    (the \code{dim}-method is sufficient).

    \item If a \code{map} is provided when \code{initialize()}ing an
    \code{"sts"} object, it is now verified that all \code{observed}
    regions are part of the \code{map} (matched by \code{row.names}).

    \item In \code{stsplot_space()}, extra (unobserved) regions of the
    \code{map} are no longer dropped but shown with a dashed border by
    default.
  }
}

}



\section{Changes in surveillance version 1.8-1 (2014-10-29)}{

\subsection{NEW FEATURES}{
  \itemize{
    \item The \code{R0}-method for \code{"twinstim"} gained an argument
    \code{newcoef} to simplify computation of reproduction numbers with
    a different parameter vector (also used for Monte Carlo CI's).

    \item New plot \code{type="neweights"} for \code{"hhh4"} fits.

    \item The \code{scores()} function allows the selection of multiple
    \code{units} (by index or name) for which to compute (averaged) proper
    scores. Furthermore, one can now select \code{which} scores to compute.

    \item Added a \code{formula}-method for \code{"hhh4"} fits to
    extract the \code{f} specifications of the three components from the
    control list.

    \item The \code{update()}-method for fitted \code{"hhh4"} models
    gained an argument \code{S} for convenient modification of component
    seasonality using \code{addSeason2formula()}.

    \item The new auxiliary function \code{layout.labels()} generates an
    \code{sp.layout} item for \code{spplot()} in order to draw labels.

    \item When generating the \code{pit()} histogram with a single
    predictive CDF \code{pdistr}, the \code{\dots} arguments can now be
    \code{x}-specific and are recycled if necessary using \code{mapply()}.
    If \code{pdistr} is a list of CDFs, \code{pit()} no longer requires
    the functions to be vectorized.

    \item New method \code{as.epidata.data.frame()}, which constructs the
    start/stop SIR event history format from a simple individual-based
    data frame (e.g., \code{hagelloch.df}).

    \item New argument \code{w} in \code{as.epidata.default()} to
    generate covariate-based weights for the force of infection in
    \code{twinSIR}. The \code{f} argument is for distance-based weights.

    \item The result of \code{profile.twinSIR()} gained a class and an
    associated \code{plot}-method.
  }
}

\subsection{MAJOR CHANGES}{
  \itemize{
    \item For multivariate \code{oneStepAhead()} predictions,
    \code{scores(..., individual=TRUE)} now returns a 3d array instead
    of a collapsed matrix. Furthermore, the scores computed by default
    are \code{c("logs","rps","dss","ses")}, excluding the normalized
    squared error score \code{"nses"} which is improper.

    \item The plot-\code{type="season"} for \code{"hhh4"} fits now by
    default plots the multiplicative effect of seasonality on the
    respective component (new argument \code{intercept=FALSE}).
    The default set of components to plot has also changed.

    \item When \code{as.epidata()} and \code{simEpidata()} calculate
    distance-based epidemic weights from the \code{f} functions, they no
    longer set the distance of an infectious individual to itself
    artificially to \code{Inf}.
    This changes the corresponding columns in the \code{"epidata"} in
    rows of currently infectious individuals, but the \code{twinSIR}
    model itself is invariant, since only rows with \code{atRiskY=1}
    contribute to the likelihood.

    \item Several modifications and corrections in \code{data("hagelloch")}.
  }
}

\subsection{MINOR CHANGES}{
  \itemize{
    \item Better plotting of \code{stsNC} objects by writing an own plot
    method for them. Prediction intervals are now shown jointly with the
    point estimate.

    \item Reduced package size by applying \code{tools::resaveRdaFiles}
    to some large datasets and by building the package with
    \code{--compact-vignettes=both}, i.e., using additional GhostScript
    compression with ebook quality, see \code{?tools::compactPDF}.

    \item Added \code{units} argument to \code{stsplot_time} to select
    only a subset of the multivariate time series for plotting.

    \item The \code{untie}-method for class \code{"epidataCS"} gained an
    argument \code{verbose} which is now \code{FALSE} by default.

    \item \code{"epidataCS"} objects store the \code{clipper} used
    during generation as attribute of \code{$events$.influenceRegion}.

    \item In \code{plotHHH4_fitted()}, the argument \code{legend.observed}
    now defaults to \code{FALSE}.

    \item The default weights for the spatio-temporal component in
    \code{hhh4} models now are \code{neighbourhood(stsObj) == 1}.
    The previous default \code{neighbourhood(stsObj)} does not make
    sense for the newly supported \code{nbOrder} neighbourhood matrices
    (shortest-path distances). The new default makes no difference for
    (old) models with binary adjacency matrices in the neighbourhood
    slot of the \code{stsObj}.

    \item The default for nonparametric weights \code{W_np()} in
    \code{hhh4()} is now to assume zero weight for neighbourhood orders
    above \code{maxlag}, i.e., \code{W_np()}'s argument \code{to0} now
    defaults to \code{TRUE}.

    \item Added a \code{verbose} argument to \code{permutationTest()},
    which defaults to \code{FALSE}. The previous behaviour corresponds
    to \code{verbose=TRUE}.

    \item \code{simulate.twinstim()} now by default uses the original
    \code{data$W} as observation region.

    \item The \code{data("measlesWeserEms")} contain two additional
    variables in the \code{@map@data} slot: \code{"vaccdoc.2004"} and
    \code{"vacc1.2004"}.

    \item The plot-method for \code{"epidata"} objects now uses colored
    lines by default.

    \item The \pkg{surveillance} package now depends on \R >= 3.0.2,
    which, effectively, is the minimum version required since
    \pkg{surveillance} 1.7-0 (see the corresponding NEWS below).

    \item The two diagnostic plots of \code{checkResidualProcess()} are
    now by default plotted side by side (\code{mfrow=c(1,2)}) instead of
    one below the other.
  }
}

\subsection{BUG FIXES}{
  \itemize{
    \item In \code{farringtonFlexible} alarms are now for
    \code{observed>upperbound} and not for \code{observed>=upperbound}
    which was not correct.

    \item Fixed duplicate \code{"functions"} element resulting from
    \code{update.twinstim(*,model=TRUE)} and ensured that
    \code{"twinstim"} objects always have the same components (some may
    be \code{NULL}).

    \item \code{animate.epidata} works again with the
    \CRANpkg{animation} package (\code{ani.options("outdir")} was
    removed in version 2.3)

    \item For \code{hhh4} models with random effects, \code{confint()}
    only worked if argument \code{parm} was specified.

    \item Computing one-sided AIC weights by simulation for \code{twinSIR}
    models with more than 2 epidemic covariates now is more robust (by
    rescaling the objective function for the quadratic programming
    solver) and twice as fast (due to code optimization).

    \item \code{simulate.twinstim(..., rmarks=NULL)} can now handle the
    case where \code{data} has no events within the simulation period
    (by sampling marks from all of \code{data$events}).

    \item The \code{lambda.h} values of simulated events in
    \code{"simEpidataCS"} objects were wrong if the model contained an
    endemic intercept (which is usually the case).

    \item Automatic choice of color breaks in the \code{animate}-method
    for class \code{"sts"} now also works for incidence maps (i.e., with
    a \code{population} argument).

    \item \code{hhh4()} did not allow the use of nonparametric
    neighbourhood weights \code{W_np()} with \code{maxlag=2}.

    \item \code{scores()} did not work for multivariate \code{oneStepAhead()}
    predictions if both \code{individual=TRUE} and \code{sign=TRUE}, and
    it could not handle a \code{oneStepAhead()} prediction of only one
    time point. Furthermore, the \code{"sign"} column of
    \code{scores(..., sign=TRUE)} was wrong (reversed).

    \item For \code{"epidataCS"} with only one event,
    \code{epidataCSplot_space()} did not draw the point.

    \item The trivial (identity) call
    \code{aggregate(stsObj, nfreq=stsObj@freq)} did not work.
  }
}

}



\section{Changes in surveillance version 1.8-0 (2014-06-16)}{

\subsection{PACKAGE INFRASTRUCTURE}{
  \itemize{
    \item Package \pkg{surveillance} now depends on newer versions of
    packages \CRANpkg{sp} (>= 1.0-15), \CRANpkg{polyCub} (>= 0.4-2),
    and \CRANpkg{spatstat} (>= 1.36-0).
    The \R packages \pkg{INLA} and \CRANpkg{runjags} are now suggested
    to support a new outbreak detection algorithm (\code{boda()}) and
    the new \code{nowcast()}ing procedure, respectively.
    The \R packages for \CRANpkg{lattice}, \CRANpkg{grid},
    \CRANpkg{gridExtra}, and \CRANpkg{scales} are suggested for
    added visualization facilities.

    \item More tests have been implemented to ensure package integrity.
    We now use \CRANpkg{testthat} instead of the outdated package
    \CRANpkg{RUnit}.

    \item \code{hhh4()} fits now have class \code{"hhh4"} instead of
    \code{"ah4"}, for consistency with \code{twinstim()},
    \code{twinSIR()}, and to follow the common convention (cp. \code{lm()}).
    Standard S3-methods for the old \code{"ah4"} name are still
    available for backwards compatibility but may be removed in the
    future.

    \item Plot variants for \code{"sts"} objects have been cleaned up:
    The functions implementing the various plot types
    (\code{stsplot_*}, previously named \code{plot.sts.*})
    are now exported and documented separately.
  }
}

\subsection{NEW FEATURES}{
  \itemize{
    \item The \code{nowcast} procedure has been completely re-written to
    handle the inherit right-truncation of reporting data (best
    visualized as a reporting triangle). The new code implements the
    generalized-Dirichlet and the hierarchical Bayesian approach described in
    \enc{Höhle}{Hoehle} and an der Heiden (2014). No backwards
    compatibility to the old nowcasting procedure is given.

    \item The package contains a new monitoring function
    \code{boda}. This is a first experimental surveillance
    implementation of the Bayesian Outbreak Detection Algorithm (BODA)
    proposed in Manitz and \enc{Höhle}{Hoehle} (2012). The function
    relies on the non-CRAN package \pkg{INLA}, which has to be installed first
    in order to use this function. Expect initial problems.


    \item New \code{toLatex}-method for \code{"sts"} objects.

    \item The new function \code{stsplot_space()} provides an improved
    map plot of disease incidence for \code{"sts"} objects aggregated
    over time. It corresponds to the new \code{type = observed ~ unit}
    of the \code{stsplot}-method, and supersedes
    \code{type = observed ~ 1|unit} (except for alarm shading).

    \item An \code{animate()}-method for the \code{"sts"} class provides
    a new implementation for animated maps (superseding the \code{plot}
    \code{type=observed ~ 1 | unit * time}) with an optional evolving
    time series plot below the map.

    \item The \code{plot()} method for \code{"sts"} objects with epochs as
    dates is now made more flexible by introducing the arguments
    \code{xaxis.tickFreq}, \code{xaxis.labelFreq} and
    \code{xaxis.labelFormat}. These allow the specification of
    tick-marks and labelling based on \code{strftime} compatible
    conversion codes -- independently if data are daily, weekly, monthly,
    etc. As a consequence, the old argument \code{xaxis.years} is
    removed. See \code{stsplot_time()} for more information.


    \item Inference for neighbourhood weights in \code{hhh4()} models:
    \code{W_powerlaw()} and \code{W_np()} both implement weights
    depending on the order of neighbourhood between regions, a power-law
    decay and nonparametric weights, i.e., unconstrained estimation of
    individual weights for each neighbourhood order.

    \item \code{hhh4()} now allows the inclusion of multiplicative
    offsets also in the epidemic components \code{"ar"} and \code{"ne"}.

    \item \code{hhh4()} now has support for \code{lag != 1} in the
    autoregressive and neighbor-driven components. The applied lags are
    stored as component \code{"lags"} of the return value (previously
    there was an unused component \code{"lag"} which was always 1 and
    has been removed now).

    \item \code{oneStepAhead()}:
    \itemize{
      \item Added support for parallel computation of predictions using
      \code{mclapply()} from package \pkg{parallel}.
      \item New argument \code{type} with a new \code{type}
      \code{"first"} to base all subsequent one-step-ahead predictions
      on a single initial fit.
      \item Nicer interpretation of \code{verbose} levels, and
      \code{txtProgressBar()}.
    }

    \item The \code{plot()}-method for fitted \code{hhh4()} objects now
    offers three additional types of plots: component seasonality,
    seasonal or time course of the dominant eigenvalue, and maps
    of estimated random intercepts. It is documented and more customizable.
    Note that argument order and some names have changed:
    \code{i} -> \code{units}, \code{title} -> \code{names}.

    \item (Deviance) \code{residuals()}-method for fitted \code{hhh4()}
    models.

    \item Added methods of \code{vcov()} and \code{nobs()}
    for the \code{"hhh4"} class. For \code{AIC()} and \code{BIC()}, the
    default methods work smoothly now (due to changes to
    \code{logLik.hhh4()} documented below).


    \item New predefined interaction functions for \code{twinstim()}:
    \code{siaf.student()} implements a \eqn{t}-kernel for the distance
    decay, and \code{siaf.step()} and \code{tiaf.step()} provide step
    function kernels (which may also be invoked by specifying the
    vector of knots as the \code{siaf} or \code{tiaf} argument in
    \code{twinstim}).

    \item Numerical integration over polygonal domains in the \code{F}
    and \code{Deriv} components of \code{siaf.powerlaw()} and
    \code{siaf.powerlawL()} is much faster and more accurate now since
    we use the new \code{polyCub.iso()} instead of \code{polyCub.SV()}
    from package \CRANpkg{polyCub}.

    \item New \code{as.stepfun()}-method for \code{"epidataCS"} objects.

    \item \code{plot.epidataCS()}:
    \itemize{
      \item The spatial plot has new arguments to automatically add
      legends to the plot: \code{legend.types} and \code{legend.counts}.
      It also gained an \code{add} argument.
      \item The temporal plot now supports type-specific sub-histograms,
      additional lines for the cumulative number of events, and an
      automatic legend.
    }

    \item The new function \code{glm_epidataCS()} can be used to fit
    an endemic-only \code{twinstim()} via \code{glm()}.
    This is mainly provided for testing purposes since wrapping into
    \code{glm} usually takes longer.
  }
}

\subsection{MAJOR CHANGES}{
  \itemize{
    \item Fitted \code{hhh4()} objects no longer contain the associated
    \code{"sts"} data twice: it is now only stored as \code{$stsObj}
    component, the hidden duplicate in \code{$control$data$.sts} was
    dropped, which makes fitted objects substantially smaller.

    \item \code{logLik.hhh4()} always returns an object of class
    \code{"logLik"} now; for random effects models, its \code{"df"}
    attribute is \code{NA_real_}. Furthermore, for non-convergent fits,
    \code{logLik.hhh4()} gives a warning and returns \code{NA_real_};
    previously, an error was thrown in this case.

    \item \code{oneStepAhead()}:
    \itemize{
      \item Default of \code{tp[2]} is now the penultimate time point of
      the fitted subset (not of the whole \code{stsObj}).
      \item \code{+1} on rownames of \code{$pred} (now the same as for
      \code{$observed}).
    }

    \item The optional \code{"twinstim"} result components
    \code{fisherinfo}, \code{tau}, and \code{functions} are always
    included. They are set to \code{NULL} if they are not applicable
    instead of missing completely (as before), such that all
    \code{"twinstim"} objects have the same list structure.

    \item \code{iafplot()} ...
    \itemize{
      \item invisibly returns a matrix containing the plotted
      values of the (scaled) interaction function (and the confidence
      interval as an attribute).
      Previously, nothing (\code{NULL}) was returned.
      \item detects a type-specific interaction function and by default
      uses \code{types=1} if it is not type-specific.
      \item has better default axis ranges.
      \item adapts to the new step function kernels (with new arguments
      \code{verticals} and \code{do.points}).
      \item supports logarithmic axes (via new \code{log} argument
      passed on to \code{plot.default}).
      \item optionally respects \code{eps.s} and \code{eps.t},
      respectively (by the new argument \code{truncated}).
      \item now uses \code{scaled=TRUE} by default.
    }

    \item The argument \code{colTypes} of
    \code{plot.epidataCS(,aggregate="space")} is deprecated (use
    \code{points.args$col} instead).

    \item The events in an \code{"epidataCS"} object no longer have
    a reserved \code{"ID"} column.
  }
}

\subsection{MINOR CHANGES}{
  \itemize{
    \item \code{hhh4()} now stores the runtime just like \code{twinstim()}.
    \item Take \code{verbose=FALSE} in \code{hhh4()} more seriously.
    \item \code{hhh4()} issues a \code{warning()} if non-convergent.
    \item The following components of a \code{hhh4()} fit now have names:
    \code{"se"}, \code{"cov"}, \code{"Sigma"}.
    \item The new default for \code{pit()} is to produce the plot.

    \item The \code{twinstim()} argument \code{cumCIF} now defaults to
    \code{FALSE}.
    \item \code{update.twinstim()} no longer uses recursive
    \code{modifyList()} for the \code{control.siaf} argument. Instead,
    the supplied new list elements (\code{"F"}, \code{"Deriv"})
    completely replace the respective elements from the original
    \code{control.siaf} specification.
    \item \code{siaf.lomax()} is now defunct (it has been deprecated
    since version 1.5-2); use \code{siaf.powerlaw()} instead.
    \item Allow the default \code{adapt}ive bandwidth to be specified via the
    \code{F.adaptive} argument in \code{siaf.gaussian()}.
    \item Unsupported options (\code{logpars=FALSE},
    \code{effRangeProb}) have been dropped from \code{siaf.powerlaw()}
    and \code{siaf.powerlawL()}.
    \item More rigorous checking of \code{tiles} in
    \code{simulate.twinstim()} and \code{intensityplot.twinstim}.

    \item \code{as.epidataCS()} gained a \code{verbose} argument.
    \item \code{animate.epidataCS()} now by default does not draw
    influence regions (\code{col.influence=NULL}), is \code{verbose} if
    \code{interactive()}, and ignores \code{sleep} on non-interactive
    devices.

    \item The \code{multiplicity}-generic and its default method have
    been integrated into \CRANpkg{spatstat} and are imported from there.
  }
}

\subsection{DATA}{
  \itemize{
    \item The polygon representation of Germany's districts (
    \code{system.file("shapes", "districtsD.RData", package="surveillance")}
    ) has been simplified further. The union of \code{districtsD} is
    used as observation window \code{W} in \code{data("imdepi")}. The
    exemplary \code{twinstim()} fit \code{data("imdepifit")} has been
    updated as well. Furthermore, \code{row.names(imdepi$events)} have
    been reset (chronological index), and numerical differences
    in \code{imdepi$events$.influenceRegion} are due to changes in
    \CRANpkg{polyclip} 1.3-0.

    \item The Campylobacteriosis data set \code{campyDE}, where absolute
    humidity is used as concurrent covariate to adjust the outbreak
    detection is added to the package to exemplify \code{boda()}.

    \item New \code{data("measlesWeserEms")} (of class \code{"sts"}),
    a corrected version of \code{data("measles.weser")} (of the old
    \code{"disProg"} class).
  }
}

\subsection{BUG FIXES}{
  \itemize{
    \item Fixed a bug in \code{LRCUSUM.runlength} where computations
    were erroneously always done under  the in-control parameter
    \code{mu0} instead of \code{mu}.

    \item Fixed a bug during alarm plots (\code{stsplot_alarm()}),
    where the use of \code{alarm.symbol} was ignored.

    \item Fixed a bug in \code{algo.glrnb} where the overdispersion
    parameter \code{alpha} from the automatically fitted \code{glm.nb}
    model (fitted by \code{estimateGLRNbHook}) was incorrectly taken as
    \code{mod[[1]]$theta} instead of \code{1/mod[[1]]$theta}. The error is
    due to a different parametrization of the negative binomial
    distribution compared to the parametrization in \enc{Höhle}{Hoehle}
    and Paul (2008).

    \item The score function of \code{hhh4()} was wrong when fitting
    endemic-only models to a \code{subset} including the first time
    point. This led to \dQuote{false convergence}.

    \item \code{twinstim()} did not work without an endemic offset if
    \code{is.null(optim.args$par)}.
  }
}

}



\section{Changes in surveillance version 1.7-0 (2013-11-19)}{

\subsection{SYNOPSIS}{
  \itemize{
    \item Package \CRANpkg{gpclib} is no longer necessary for the
    construction of \code{"epidataCS"}-objects. Instead, we make use of
    the new dedicated package \CRANpkg{polyclip} (licensed under the
    BSL) for polygon clipping operations (via
    \code{spatstat::intersect.owin()}). This results in a slightly
    different \code{$events$.influenceRegion} component of
    \code{"epidataCS"} objects, one reason being that
    \pkg{polyclip} uses integer arithmetic.
    Change of \code{twinstim()} estimates for a newly created
    \code{"epidataCS"} compared with the same data prepared in earlier
    versions should be very small (e.g., for \code{data("imdepifit")}
    the mean relative difference of coefficients is 3.7e-08, while the
    \code{logLik()} is \code{all.equal()}).
    As an alternative, \pkg{rgeos} can still be chosen to do the polygon
    operations.

    \item The \pkg{surveillance}-internal code now depends on
    \R >= 2.15.2 (for \code{nlminb()} \code{NA} fix of PR#15052,
    consistent \code{rownames(model.matrix)} of PR#14992,
    \code{paste0()}, \code{parallel::mcmapply()}). However, the
    required recent version of \pkg{spatstat} (1.34-0, for
    \pkg{polyclip}) actually needs \R >= 3.0.2, which therefore also
    applies to \pkg{surveillance}.

    \item Some minor new features and changes are documented below.
  }
}

\subsection{NEW FEATURES}{
  \itemize{
    \item Functions \code{unionSpatialPolygons()} and
    \code{intersectPolyCircle()} are now exported. Both are wrappers
    around functionality from different packages supporting polygon
    operations: for determining the union of all subpolygons of a
    \code{"SpatialPolygons"} object, and the intersection of a polygonal
    and a circular domain, respectively.
    \item \code{discpoly()} moved back from \CRANpkg{polyCub}
    to \pkg{surveillance}.
  }
}

\subsection{MINOR CHANGES}{
  \itemize{
    \item \pkg{surveillance} now Depends on \CRANpkg{polyCub} (>= 0.4-0)
    and not only Imports it (which avoids \code{::}-references in
    .GlobalEnv-made functions).
    \item Nicer default axis labels for \code{iafplot()}.
    \item For \code{twinstim()}, the default is now to \code{trace}
    every iteration instead of every fifth only.
    \item Slightly changed default arguments for \code{plot.epidata()}:
    \code{lwd} (1->2), \code{rug.opts} (\code{col} is set according to
    \code{which.rug})
    \item \code{twinstim()} saves the vector of \code{fixed}
    coefficients as part of the returned \code{optim.args} component,
    such that these will again be held fixed upon \code{update()}.
    \item The \code{plot}-method for \code{hhh4()}-fits allows for
    region selection by name.
  }
}

}



\section{Changes in surveillance version 1.6-0 (2013-09-03)}{

\subsection{SYNOPSIS}{
  \itemize{
    \item The \code{polyCub}-methods for cubature over polygonal domains
    have been moved to the new dedicated package \CRANpkg{polyCub},
    since they are of a rather general use. The \code{discpoly()} function
    has also been moved to that package.

    \item As a replacement for the license-restricted \pkg{gpclib} package,
    the \pkg{rgeos} package is now used by default
    (\code{surveillance.options(gpclib=FALSE)}) in generating
    \code{"epidataCS"} (polygon intersections, slightly slower).
    Therefore, when installing \pkg{surveillance} version 1.6-0, the
    system requirements for \CRANpkg{rgeos} have to be met, i.e., GEOS
    must be available on the system. On Linux variants this means
    installing \file{libgeos} (\file{libgeos-dev}).

    \item The improved Farrington method described in Noufaily et al.
    (2012) is now available as function \code{farringtonFlexible()}.

    \item New handling of reference dates in \code{algo.farrington()} for
    \code{"sts"} objects with \code{epochAsDate=TRUE}. Instead of always
    going back in time to the next Date in the \code{"epoch"} slot, the
    function now determines the \emph{closest} Date. Note that this
    might lead to slightly different results for the upperbound
    compared to previously. Furthermore, the functionality is only
    tested for weekly data (monthly data are experimental). The same
    functionality applies to \code{farringtonFlexible()}.

    \item To make the different retrospective modelling frameworks of
    the \pkg{surveillance} package jointly applicable, it is now possible
    to convert (aggregate) \code{"epidataCS"}
    (continuous-time continuous-space data) into an \code{"sts"} object
    (multivariate time series of counts) by the new function
    \code{epidataCS2sts}.

    \item Simulation from \code{hhh4} models has
    been re-implemented, which fixes a bug and makes it more flexible
    and compatible with a wider class of models.

    \item The \code{map}-slot of the \code{"sts"} class now requires
    \code{"SpatialPolygons"} (only) instead of
    \code{"SpatialPolygonsDataFrame"}.

    \item Re-implementation of \code{oneStepAhead()} for
    \code{hhh4}-models with a bug fix, some speed-up and more options.

    \item Slight speed-up for \code{hhh4()} fits,
    e.g., by more use of \code{.rowSums()} and \code{.colSums()}.

    \item Crucial speed-up for \code{twinstim()} fits by more efficient
    code: \code{mapply}, dropped clumsy \code{for}-loop in
    \code{fisherinfo}, new argument \code{cores} for parallel
    computing via forking (not available on Windows).

    \item Some further new features, minor changes, and bug fixes are
    described in the following subsections.
  }
}

\subsection{NEW FEATURES}{
  \itemize{
    \item Using \code{tiaf.exponential()} in a \code{twinstim()} now works
    with \code{nTypes=1} for multi-type data.

    \item A legend can be added automatically in \code{iafplot()}.

    \item The \code{untie} methods are now able to produce jittered points
    with a required minimum separation (\code{minsep}).

    \item \code{simulate.ah4} gained a \code{simplify} argument.

    \item New \code{update}-method for fitted \code{hhh4}-models
    (class \code{"ah4"}).

    \item \code{oneStepAhead()} has more options: specify time range
    (not only start), choose type of start values, \code{verbose}
    argument.

    \item \code{pit()} allows for a list of predictive distributions
    (\code{pdistr}), one for each observation \code{x}.

    \item New spatial auxiliary function \code{polyAtBorder()}
    indicating polygons at the border (for a \code{"SpatialPolygons"}
    object).

    \item \code{animate.epidataCS()} allows for a \code{main} title and
    can show a progress bar.
  }
}

\subsection{MINOR CHANGES}{
  \itemize{
    \item Changed parametrization of \code{zetaweights()} and completed
    its documentation (now no longer marked as experimental).

    \item \code{twinstim(...)$converged} is \code{TRUE} if
    the optimization routine converged (as before) but contains
    the failure message otherwise.

    \item Increased default \code{maxit} for the Nelder-Mead optimizer
    in \code{hhh4} from 50 to 300, and removed default artificial lower
    bound (-20) on intercepts of epidemic components.

    \item Renamed returned list from \code{oneStepAhead}
    (mean->pred, x->observed, params->coefficients,
    variances->Sigma.orig) for consistency, and
    \code{oneStepAhead()$psi} is only non-\code{NULL} if we have a
    NegBin model.

    \item Argument order of \code{pit()} has changed, which is also
    faster now and got additional arguments \code{relative} and
    \code{plot}.

    \item \code{twinstim(...)$runtime} now contains the complete
    information from \code{proc.time()}.
  }
}

\subsection{BUG FIXES}{
  \itemize{
    \item Fixed a bug in function
    \code{refvalIdxByDate()} which produced empty reference values
    (i.e. \code{NA}s) in case the Date entries of \code{epoch} were not
    mondays. Note: The function works by subtracting \code{1:b} years from the
    date of the range value and then takes the span \code{-w:w} around this
    value. For each value in this set it is determined whether the
    closest date in the epoch slot is obtained by going forward or
    backward. Note that this behaviour is now slightly changed compared
    to previously, where we \emph{always} went back in time.

    \item \code{algo.farrington()}: Reference values too far back in time
    and hence not being in the \code{"epoch"} slot of the \code{"sts"}
    object are now ignored (previously the resulting \code{NA}s caused the
    function to halt). A warning is displayed in this case.

    \item \code{hhh4}: The entry \eqn{(5,6)} of the marginal Fisher
    information matrix in models with random intercepts in all three
    components was incorrect.
    If \code{nlminb} was used as optimizer for the variance parameters
    (using the negative marginal Fisher information as Hessian), this
    could have caused false convergence (with warning) or minimally
    biased convergence (without warning).
    As a consequence, the \code{"Sigma.cov"} component of the
    \code{hhh4()} result, which is the inverse of the marginal Fisher
    information matrix at the MLE, was also wrong.

    \item \code{untie.matrix()} could have produced jittering greater than
    the specified \code{amount}.

    \item \code{hhh4}: if there are no random intercepts, the
    redundant \code{updateVariance} steps are no longer evaluated.

    \item \code{update.twinstim()} did not work with
    \code{optim.args=..1} (e.g., if updating a list of models with lapply).
    Furthermore, if adding the \code{model} component only, the
    \code{control.siaf} and \code{optim.args} components were lost.

    \item \code{earsC} should now also work with multivariate
    \code{sts} time-series objects.

    \item The last week in \code{data(fluBYBW)} (row index 417) has been
    removed. It corresponded to week 1 in year 2009 and was wrong
    (an artifact, filled with zero counts only).
    Furthermore, the regions in \code{@map} are now ordered the same as
    in \code{@observed}.

    \item Fixed start value of the overdispersion parameter in
    \code{oneStepAhead} (must be on internal log-scale, not
    reparametrized as returned by \code{coef()} by default).

    \item When subsetting \code{"sts"} objects in time, \code{@start} was
    updated but not \code{@epoch}.

    \item \code{pit} gave \code{NA} results if any \code{x[-1]==0}.

    \item The returned \code{optim.args$par} vector in \code{twinstim()}
    was missing any fixed parameters.

    \item \code{hhh4()} did not work with time-varying neighbourhood
    weights due to an error in the internal \code{checkWeightsArray()}
    function.
  }
}

}



\section{Changes in surveillance version 1.5-4 (2013-04-21)}{

\subsection{SYNOPSIS}{
  \itemize{
    \item Fixed obsolete \code{.path.package()} calls.
    \item Small corrections in the documentation.
    \item \code{update.twinstim()} performs better in preserving
    the original initial values of the parameters.
    \item New pre-defined spatial interaction function
    \code{siaf.powerlawL()}, which implements a _L_agged power-law
    kernel, i.e. accounts for uniform short-range dispersal.
  }
}

}



\section{Changes in surveillance version 1.5-2 (2013-03-15)}{

\subsection{SYNOPSIS}{
  \itemize{
    \item New method for outbreak detection: \code{earsC}
    (CUSUM-method described in the CDC Early Aberration Reporting
    System, see Hutwagner et al, 2003).
    \item New features and minor bug fixes for the "\code{twinstim}"
    part of the package (see below).
    \item Yet another p-value formatting function \code{formatPval()}
    is now also part of the \pkg{surveillance} package.
    \item \code{polyCub.SV()} now also accepts objects of classes
    \code{"Polygon"} and \code{"Polygons"} for convenience.
    \item \code{siaf.lomax} is deprecated and replaced by
    \code{siaf.powerlaw} (re-parametrization).
  }
}

\subsection{NEW FEATURES (\code{twinstim()}-related)}{
  \itemize{
    \item The temporal \code{plot}-method for class \code{"epidataCS"}
    now understands the \code{add} parameter to add the histogram to an
    existing plot window, and auto-transforms the \code{t0.Date}
    argument using \code{as.Date()} if necessary.
    \item \code{nobs()} methods for classes \code{"epidataCS"} and
    \code{"twinstim"}.
    \item New argument \code{verbose} for \code{twinstim()} which, if
    set to \code{FALSE}, disables the printing of information messages
    during execution.
    \item New argument \code{start} for \code{twinstim()}, where (some)
    initial parameter values may be provided, which overwrite those in
    \code{optim.args$par}, which is no longer required (as a naive
    default, a crude estimate for the endemic intercept and zeroes for
    the other parameters are used).
    \item Implemented a wrapper \code{stepComponent()} for \code{step()}
    to perform algorithmic component-specific model selection in
    \code{"twinstim"} models. This also required the implementation of
    suitable \code{terms()} and \code{extractAIC()} methods. The single-step
    methods \code{add1()} and \code{drop1()} are also available.
    \item The \code{update.twinstim()} method now by default uses the
    parameter estimates from the previous model as initial values for
    the new fit (new argument \code{use.estimates = TRUE}).
    \item \code{as.epidataCS()} checks for consistency of the area of
    \code{W} and the (now really obligatory) area column in
    \code{stgrid}.
    \item \code{simulate.twinstim()} now by default uses the previous
    \code{nCircle2Poly} from the \code{data} argument.
    \item \code{direction} argument for \code{untie.epidataCS()}.
    \item The \code{toLatex}-method for \code{"summary.twinstim"} got
    different defaults and a new argument \code{eps.Pvalue}.
    \item New \code{xtable}-method for \code{"summary.twinstim"} for
    printing the covariate effects as risk ratios (with CI's and p-values).
  }
}

\subsection{NEW FEATURES (\code{hhh4()}-related)}{
  \itemize{
    \item New argument \code{hide0s} in the \code{plot}-method for class
    \code{"ah4"}.
    \item New argument \code{timevar} for \code{addSeason2formula()},
    which now also works for long formulae.
  }
}

}



\section{Changes in surveillance version 1.5-1 (2012-12-14)}{

\subsection{SYNOPSIS}{
  \itemize{
    \item The \pkg{surveillance} package is again backward-compatible
    with \R version 2.14.0, which is now declared as the minimum
    required version.
  }
}

}



\section{Changes in surveillance version 1.5-0 (2012-12-12)}{

\subsection{SYNOPSIS}{
  \itemize{
    \item This new version mainly improves upon the \code{twinstim()} and
    \code{hhh4()} implementations (see below).
    \item As requested by the CRAN team, examples now run faster. Some
    are conditioned on the value of the new package option
    \code{"allExamples"}, which usually defaults to \code{TRUE} (but is
    set to \code{FALSE} for CRAN checking, if timings are active).
    \item Moved some rarely used package dependencies to
    \dQuote{Suggests:}, and also removed some unused packages from there.
    \item Dropped strict dependence on \CRANpkg{gpclib}, which has a
    restricted license, for the \pkg{surveillance} package to be clearly
    GPL-2. Generation of \code{"epidataCS"} objects, which makes use of
    \pkg{gpclib}'s polygon intersection capabilities, now requires prior
    explicit acceptance of the \pkg{gpclib} license via setting
    \code{surveillance.options(gpclib = TRUE)}. Otherwise,
    \code{as.epidataCS()} and \code{simEpidataCS()} may not be used.
  }
}

\subsection{NEW FEATURES (\code{twinstim()}-related)}{
  \itemize{
    \item Speed-up by memoisation of the \code{siaf} cubature (using
    the \CRANpkg{memoise} package).
    \item Allow for \code{nlm}-optimizer (really not recommended).
    \item Allow for \code{nlminb}-specific control arguments.
    \item Use of the expected Fisher information matrix can be disabled
    for \code{nlminb} optimization.
    \item Use of the \code{effRange}-trick can be disabled in
    \code{siaf.gaussian()} and \code{siaf.lomax()}. The default
    \code{effRangeProb} argument for the latter has been changed from
    0.99 to 0.999.
    \item The \code{twinstim()} argument \code{nCub} has been replaced
    by the new \code{control.siaf} argument list. The old
    \code{nCub.adaptive} indicator became a feature of the
    \code{siaf.gaussian()} generator (named \code{F.adaptive} there) and
    does no longer depend on the \code{effRange} specification, but uses
    the bandwidth \code{adapt*sd}, where the \code{adapt} parameter may be
    specified in the \code{control.siaf} list in the \code{twinstim()}
    call. Accordingly, the components \code{"nCub"} and
    \code{"nCub.adaptive"} have been removed from the result
    of \code{twinstim()}, and are replaced by \code{"control.siaf"}.
    \item The \code{"method"} component of the \code{twinstim()} result
    has been replaced by the whole \code{"optim.args"}.
    \item The new \code{"Deriv"} component of \code{siaf} specifications
    integrates the \dQuote{siaf$deriv} function over a polygonal domain.
    \code{siaf.gaussian()} and \code{siaf.lomax()} use \code{polyCub.SV()}
    (with intelligent \code{alpha} parameters) for this task
    (previously: midpoint-rule with naive bandwidth)
    \item \code{scaled} \code{iafplot()} (default \code{FALSE}). The
    \code{ngrid} parameter has been renamed to \code{xgrid} and is more
    general.
    \item The \code{"simulate"} component of \code{siaf}'s takes an
    argument \code{ub} (upperbound for distance from the source).
    \item Numerical integration of spatial interaction functions with an
    \code{Fcircle} trick is more precise now; this slightly changes
    previous results.
    \item New \acronym{S3}-generic \code{untie()} with a method for the
    \code{"epidataCS"} class (to randomly break tied event times and/or
    locations).
    \item Renamed \code{N} argument of \code{polyCub.SV()} to
    \code{nGQ}, and \code{a} to \code{alpha}, which both have new
    default values.
    The optional polygon rotation proposed by Sommariva &
    Vianello is now also implemented (based on the corresponding MATLAB
    code) and available as the new \code{rotation} argument.
    \item The \code{scale.poly()} method for \code{"gpc.poly"} is now
    available as \code{scale.gpc.poly()}. The default return class of
    \code{discpoly()} was changed from \code{"gpc.poly"} to
    \code{"Polygon"}.
    \item An \code{intensityplot()}-method is now also implemented for
    \code{"simEpidataCS"}.
  }
}

\subsection{NEW FEATURES (\code{hhh4()}-related)}{
  \itemize{
    \item Significant speed-up (runs about 6 times faster now, amongst
    others by many code optimizations and by using sparse \CRANpkg{Matrix}
    operations).
    \item \code{hhh4()} optimization routines can now be customized for
    the updates of regression and variance parameters seperately, which
    for instance enables the use of Nelder-Mead for the variance
    updates, which seems to be more stable/robust as it does
    not depend on the inverse Fisher info and is usually faster.
    \item The \code{ranef()} extraction function for \code{"ah4"} objects
    gained a useful \code{tomatrix} argument, which re-arranges random
    effects in a unit x effect matrix (also transforming CAR effects
    appropriately).
    \item Generalized \code{hhh4()} to also capture parametric
    neighbourhood weights (like a power-law decay). The new function
    \code{nbOrder()} determines the neighbourhood order matrix
    from a binary adjacency matrix (depends on package \CRANpkg{spdep}).
    \item New argument \code{check.analyticals} (default \code{FALSE})
    mainly for development purposes.
  }
}

\subsection{BUG FIXES}{
  \itemize{
    \item Fixed sign of observed Fisher information matrix in
    \code{twinstim}.
    \item Simulation from the Lomax kernel is now correct (via polar
    coordinates).
    \item Fixed \code{modifyListcall()} to also work with updated \code{NULL} arguments.
    \item Fixed wrong Fisher information entry for the overdispersion
    parameter in \code{hhh4}-models.
    \item Fixed wrong entries in penalized Fisher information wrt the
    combination fixed effects x CAR intercept.
    \item Fixed indexing bug in penalized Fisher calculation in the case
    of multiple overdispersion parameters and random intercepts.
    \item Fixed bug in Fisher matrix calculation concerning the relation
    of unit-specific and random effects (did not work previously).
    \item Improved handling of non-convergent / degenerate solutions during
    \code{hhh4} optimization. This involves using \code{ginv()} from
    package \CRANpkg{MASS}, if the penalized Fisher info is singular.
    \item Correct labeling of overdispersion parameter in
    \code{"ah4"}-objects.
    \item Some control arguments of \code{hhh4()} have more clear
    defaults.
    \item The result of \code{algo.farrington.fitGLM.fast()} now
    additionally inherits from the \code{"lm"} class to avoid warnings
    from \code{predict.lm()} about fake object.
    \item Improved \file{NAMESPACE} imports.
    \item Some additional tiny bug fixes, see the subversion log on
    R-Forge for details.
  }
}

}



\section{Changes in surveillance version 1.4-2 (2012-08-17)}{

  \subsection{SYNOPSIS}{
    \itemize{
      \item This is mainly a patch release for the \code{twinstim}-related
      functionality of the package.
      \item Apart from that, the package is now again compatible with older
      releases of \R (< 2.15.0) as intended (by defining \code{paste0()} in
      the package namespace if it is not found in \R \pkg{base} at
      installation of the \pkg{surveillance} package).
    }
  }

  \subsection{NEW FEATURES}{
    \itemize{
      \item Important new \code{twinstim()}-feature: fix parameters
      during optimization.
      \item Useful \code{update}-method for \code{"twinstim"}-objects.
      \item New \code{[[}- and \code{plot}-methods for
      \code{"simEpidataCSlist"}-objects.
      \item \code{simEpidataCS()} received tiny bug fixes and is now
      able to simulate from epidemic-only models.
      \item \code{R0}-method for \code{"simEpidataCS"}-objects (actually
      a wrapper for \code{R0.twinstim()}).
      \item Removed \code{dimyx} and \code{eps} arguments from
      \code{R0.twinstim()}; now uses \code{nCub} and
      \code{nCub.adaptive} from the fitted model and applies the same
      (numerical) integration method.
      \item \code{animate.epidata} is now compatible with the
      \CRANpkg{animation} package.
      \item More thorough documentation of \code{"twinstim"}-related
      functions \emph{including many examples}.
    }
  }

  \subsection{BUG FIXES (\code{"twinstim"}-related)}{
    \itemize{
      \item \code{nlminb} (instead of \code{optim}'s \code{"BFGS"}) is
      now the default optimizer (as already documented).
      \item The \code{twinstim}-argument \code{nCub} can now be omitted when
      using \code{siaf.constant()} (as documented) and is internally set to
      \code{NA_real_} in this case. Furthermore, \code{nCub} and
      \code{nCub.adaptive} are set to \code{NULL} if there is
      no epidemic component in the model.
      \item \code{toLatex.summary.twinstim} now again works for
      \code{summary(*, test.iaf=FALSE)}.
      \item \code{print}- and \code{summary}-methods for
      \code{"epidataCS"} no longer assume that the \code{BLOCK} index
      starts at 1, which may not be the case when using a subset in
      \code{simulate.twinstim()}.
      \item The \code{"counter"} step function returned by
      \code{summary.epidataCS()} does no longer produce false
      numbers of infectives (they were lagged by one timepoint).
      \item \code{plot.epidataCS()} now resolves \dots correctly and
      the argument \code{colTypes} takes care of a possible
      \code{subset}.
      \item \code{simEpidataCS()} now also works for endemic-only models
      and is synchronised with \code{twinstim()} regarding the
      way how \code{siaf} is numerically integrated (including the
      argument \code{nCub.adaptive}).
      \item Fixed problem with \code{simEpidataCS()} related to missing
      \file{NAMESPACE} imports (and re-exports) of \code{marks.ppp} and
      \code{markformat.default} from \CRANpkg{spatstat}, which are required
      for \code{spatstat::runifpoint()} to work, probably because
      \pkg{spatstat} currently does not register its S3-methods.
      \item Improved error handling in \code{simEpidataCS()}. Removed a
      \code{browser()}-call and avoid potentially infinite loop.
    }
  }

  \subsection{BUG FIXES (\code{"twinSIR"}-related)}{
    \itemize{
      \item The \code{.allocate} argument of \code{simEpidata()} has
      now a fail-save default.
      \item Simulation without endemic \code{cox()}-terms now works.
    }
  }

  \subsection{MINOR CHANGES}{
    \itemize{
      \item Simplified \code{imdepi} data to monthly instead of weekly
      intervals in \code{stgrid} for faster examples and reduced package
      size.
      \item The environment of all predefined interaction functions for
      \code{twinstim()} is now set to the \code{.GlobalEnv}. The previous
      behaviour of defining them in the \code{parent.frame()} could have
      led to huge \code{save()}'s of \code{"twinstim"} objects even with
      \code{model=FALSE}.
      \item \code{simulate.twinSIR} only returns a list of epidemics if
      \code{nsim > 1}.
      \item \code{simulate.twinstim} uses \code{nCub} and
      \code{nCub.adaptive} from fitted object as defaults.
      \item Removed the \dots-argument from \code{simEpidataCS()}.
      \item The coefficients returned by \code{simEpidataCS()} are now stored
      in a vector rather than a list for compatibility with
      \code{"twinstim"}-methods.
      \item Argument \code{cex.fun} of \code{intensityplot.twinstim()} now
      defaults to the \code{sqrt} function (as in \code{plot.epidataCS()}.
    }
  }
}


\section{Changes in surveillance version 1.4 (2012-07-26)}{

  \subsection{SYNOPSIS}{
    \itemize{
      \item Besides minor bug fixes, additional functionality has entered the package
            and a new attempt is made to finally release a new version on CRAN
            (version 1.3 has not appeared on CRAN), including a proper \file{NAMESPACE}.
    }
  }

  \subsection{NEW FEATURES}{
    \itemize{
      \item Support for non-parametric back-projection using the function
            \code{backprojNP()} which returns an object of the new
	    \code{"stsBP"} class which inherits from \code{"sts"}.
      \item Bayesian nowcasting for discrete time count data is implemented in
            the function \code{nowcast()}.
      \item Methods for cubature over polygonal domains can now also visualize what
            they do.  There is also a new quasi-exact method for cubature of the
            bivariate normal density over polygonal domains.  The
	    function \code{polyCub()} is a wrapper for the different
	    methods.
      \item \code{residuals.twinstim()} and \code{residuals.twinSIR()}:
	    extract the \dQuote{residual process}, see \cite{Ogata
            (1988)}. The residuals of \code{"twinSIR"} and
            \code{"twinstim"} models may be checked graphically by the new
            function \code{checkResidualProcess()}.
      \item Many new features for the \code{"twinstim"} class of
            self-exciting spatio-temporal point process models (see
	    below).
    }
  }

  \subsection{NEW FEATURES AND SIGNIFICANT CHANGES FOR \code{"twinstim"}}{
    \itemize{
      \item Modified arguments of \code{twinstim()}: new ordering, new
            argument \code{nCub.adaptive}, removed argument
	    \code{typeSpecificEndemicIntercept} (which is now specified as part of
	    the \code{endemic} formula as \code{(1|type)}).
      \item Completely rewrote the \code{R0}-method (calculate \dQuote{trimmed} and
            \dQuote{untrimmed} \eqn{R_0} values)
      \item The \dQuote{trimmed} \code{R0} values are now part of the
            result of the model fit, as well as \code{bbox(W)}. The
	    model evaluation environment is now set as attribute of the
	    result if \code{model=TRUE}.
      \item New predefined spatial kernel: the Lomax power law kernel
            \code{siaf.lomax()}
      \item \code{plot}-methods for \code{"twinstim"}
            (\code{intensityplot()} and \code{iafplot()})
      \item \code{as.epidataCS()} now auto-generates the stop-column if this is missing
      \item \code{print}-method for class \code{"summary.epidataCS"}
      \item \code{[}- and subset-method for \code{"epidataCS"}
            (subsetting \code{...$events})
      \item \code{plot}-method for \code{"epidataCS"}
    }
  }

  \subsection{MINOR CHANGES}{
    \itemize{
      \item Improved documentation for the new functionalities.
      \item Updated references.
      \item \code{twinSIR}'s \code{intensityPlot} is now a method of the
            new S3-generic function \code{intensityplot}.
    }
  }
}


\section{Changes in surveillance version 1.3 (2011-04-25)}{
  \subsection{SYNOPSIS}{
    \itemize{
      \item This is a major realease integrating plenty of new code (unfortunately
            not all documented as good as it could be). This includes code
            for the \code{"twinstim"} and the \code{"hhh4"} model.
	    The \code{"twinSIR"} class of models has been
            migrated from package \pkg{RLadyBug} to \pkg{surveillance}.
            It may take a while before this version will become available from CRAN.
            For further details see below.
    }
  }
  \subsection{SIGNIFICANT CHANGES}{
    \itemize{
      \item Renamed the \code{"week"} slot of the \code{"sts"} S4 class to \code{"epoch"}.
            All saved data objects have accordingly be renamed, but some hazzle
            is to be expected if one you have old \code{"sts"} objects stored in binary
            form. The function \code{convertSTS()} can be used to
	    convert such \dQuote{old school} \code{"sts"} objects.
      \item Removed the functions \code{algo.cdc()} and \code{algo.rki()}.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Support for \code{"twinSIR"} models (with associated
            \code{"epidata"} objects) as described
            in \cite{Höhle (2009)} has been moved from package
	    \pkg{RLadyBug} to \pkg{surveillance}.
            That means continuous-time discrete-space \acronym{SIR} models.
      \item Support for \code{"twinstim"} models as described in
            \cite{Meyer et al (2012)}. That means continuous-time
            continuous-space infectious disease models.
      \item Added functionality for non-parametric back projection
            (\code{backprojNP()}) and
            now-casting (\code{nowcast()}) based on \code{"sts"} objects.
    }
  }
}

\section{Changes in surveillance version 1.2-2}{
  \itemize{
    \item Replaced the deprecated getSpPPolygonsLabptSlots method with calls
      to the coordinates method when plotting the map slot.
    \item Minor proof-reading of the documentation.
    \item Added an argument \code{"extraMSMargs"} to the algo.hmm function.
    \item Fixed bug in \code{outbreakP()} when having observations equal to zero
      in the beginning. Here, \eqn{\hat{\mu}^{C1}} in (5) of \cite{Frisen et al (2008)}
      is zero and hence the log-based summation in the code failed.
      Changed to product as in the original code, which however might be
      less numerically stable.
    \item Fixed bug in stcd which added one to the calculated index of idxFA and idxCC.
      Thanks to Thais Rotsen Correa for pointing this out.
  }
}

\section{Changes in surveillance version 1.2-1 (2010-06-10)}{
  \itemize{
    \item Added \code{algo.outbreakP()} (\cite{Frisen & Andersson, 2009}) providing a
      semiparametric approach for outbreak detection for Poisson
      distributed variables.
    \item Added a pure \R function for extracting ISO week and year from Date
      objects. This function (isoWeekYear) is only called if "\%G" and "\%V"
      format strings are used on Windows (\code{sessionInfo()[[1]]$os == "mingw32"})
      as this is not implemented for \code{"format.Date"} on Windows.
      Thanks to Ashley Ford,  University of Warwick, UK for identifying
      this Windows specific bug.
    \item For \code{algo.farrington()} a faster fit routine \code{"algo.farrington.fitGLM.fast"}
      has been provided by Mikko Virtanen, National Institute for Health
      and Welfare, Finland. The new function calls \code{glm.fit()}
      directly, which gives a doubling of speed for long series. However, if one
      wants to process the fitted model output some of the GLM routines might
      not work on this output. For backwards compability the argument
      \code{control$fitFun = "algo.farrington.fitGLM"} provides the old (and slow)
      behaviour.
  }
}

\section{Changes in surveillance version 1.1-6 (2010-05-25)}{
  \itemize{
    \item A few minor bug fixes
    \item Small improvements in the C-implementation of the \code{twins()}
      function by Daniel Sabanés Bové fixing the segmentation fault
      issue on 64-bit architectures.
  }
}

\section{Changes in surveillance version 1.1-2 (2009-10-15)}{
  \itemize{
    \item Added the functions categoricalCUSUM and LRCUSUM.runlength
      for the CUSUM monitoring of general categorical time series
      (binomial, beta-binomial, multinomial, ordered response,
      Bradley-Terry models).
    \item Added the functions pairedbinCUSUM and pairedbinCUSUM.runlength
      implementing the CUSUM monitoring and run-length computations for
      a paired binary outcome as described in Steiner et al. (1999).
    \item Experimental implementation of the prospective space-time cluster
      detection described in Assuncao and Correa (2009).
    \item Added a \code{demo("biosurvbook")} containing the code of an upcoming
      book chapter on how to use the surveillance package. This
      contains the description of ISO date use, negative binomial CUSUM,
      run-length computation, etc. From an applicational point of view
      the methods are illustrated by Danish mortality monitoring.
    \item Fixed a small bug in algo.cdc found by Marian Talbert Allen
      which resulted in the control$m argument being ignored.
    \item The constructor of the sts class now uses the argument
      \code{"epoch"} instead of weeks to make clearer that also
      daily, monthly or other data can be handled.
    \item Added additional epochAsDate slot to sts class. Modified
      plot functions so they can handle ISO weeks.
    \item algo.farrington now also computes quantile and median of
      the predictive distribution. Furthermore has the computation
      of reference values been modified so its a) a little bit faster
      and b) it is also able to handle ISO weeks now. The reference values
      for date t0 are calculated as follows:
      For i, i=1,..., b look at date t0 - i*year. From this date on move
      w months/weeks/days to the left and right. In case of weeks:
       For each of these
      determined time points go back in time to the closest Monday
    \item Renamed the functions obsinyear to epochInYear, which now also
      handles objects of class Date.
  }
}

\section{Changes in surveillance version 1.0-2 (2009-03-06)}{
  \itemize{
    \item Negative Binomial CUSUM or the more general NegBin likelihood ratio
      detector is now implemented as part of algo.glrnb.
      This includes the back calculation of the required number of cases
      before an alarm.
    \item Time varying proportion binomial CUSUM.
  }
}

\section{Changes in surveillance version 0.9-10}{
  \itemize{
    \item Current status: Development version available from
      \url{http://surveillance.r-forge.r-project.org/}
    \item Rewriting of the plot.sts.time.one function to use polygons
      instead of lines for the number of observed cases. Due cause
      a number of problems were fixed in the plotting of the legend.
      Plotting routine now also handles binomial data, where the
      number of observed cases y are stored in \code{"observed"} and the
      denominator data n are stored in \code{"populationFrac"}.
    \item Problems with the aggregate function not operating correctly
      for the populationFrac were fixed.
    \item The \code{"rogerson"} wrapper function for algo.rogerson was modified so it
      now works better for distribution \code{"binomial"}. Thus a time varying
      binomial cusum can be run by calling
           \code{rogerson( x, control(..., distribution="binomial"))}
    \item An experimental implementation of the twins model documented in
      Held, L., Hofmann, M., Höhle, M. and Schmid V. (2006) A two-component
      model for counts of infectious diseases, Biostatistics, 7, pp.
      422--437 is now available as algo.twins.
  }
}

\section{Changes in surveillance version 0.9-9 (2008-01-21)}{
  \itemize{
    \item Fixed a few small problems
      which gave warnings in the CRAN distribution
  }
}

\section{Changes in surveillance version 0.9-8 (2008-01-19)}{
  \itemize{
    \item The algo_glrpois function now has an additional \code{"ret"} arguments,
      where one specifies the return type. The arguments of the underlying
      c functions have been changed to include an additional direction and
      return type value arguments.
    \item added restart argument to the algo.glrpois control object, which
      allows the user to control what happens after the first alarm has been
      generated
    \item experimental algo.glrnb function is added to the package. All calls to
      algo.glrpois are now just alpha=0 calls to this function. However,
      the underlying C functions differentiate between poisson and negative case
  }
}
