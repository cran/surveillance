\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
%% some pre-defined commands: \R, \code, \acronym, \url, \file, \pkg

\name{NEWS}
\title{News for Package 'surveillance'}
\encoding{latin1}



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

    \item Simulation from \code{hhh4} models (class \code{"ah4"}) has
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
            migrated from package \CRANpkg{RLadyBug} to \pkg{surveillance}.
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
	    \CRANpkg{RLadyBug} to \pkg{surveillance}.
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
