\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
%% some pre-defined commands: \R, \code, \acronym, \url, \file, \pkg

\name{NEWS}
\title{News for Package 'surveillance'}
\encoding{latin1}

\section{Changes in surveillance version 1.4}{
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
            continuous-space infectious models .
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
