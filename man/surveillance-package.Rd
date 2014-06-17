% Note @ R 3.0.2: Unfortunately, setting \RdOpts{stage=build} globally
% does not work as expected, so we have to specify it in each \Sexpr

\encoding{latin1}
\name{surveillance-package}
\alias{surveillance-package}
\alias{surveillance}
\docType{package}

\title{
  \Sexpr[stage=build]{(meta <- packageDescription("surveillance", encoding="latin1"))$Title}
  % If !is.na(encoding), do iconv() from "latin1" (DESCRIPTION: Encoding) to 'encoding'
  % Do we really need re-encoding, i.e., encoding="latin1", to make the \Sexpr's work on Windows?
}

\description{
  A package implementing statistical methods for the modelling and
  change-point detection in time series of counts, proportions and
  categorical data, as well as for the modeling of continuous-time
  epidemic phenomena, e.g. discrete-space setups such as the spatially
  enriched Susceptible-Exposed-Infectious-Recovered (SEIR) models for
  surveillance data, or continuous-space point process data such as the
  occurrence of disease or earthquakes. Main focus is on outbreak
  detection in count data time series originating from public health
  surveillance of infectious diseases, but applications could just as well
  originate from environmetrics, reliability engineering, econometrics or
  social sciences.
  Currently the package contains implementations of typical outbreak
  detection procedures such as Stroup et. al (1989), Farrington et al
  (1996), Rossi et al (1999), Rogerson and Yamada (2001), a Bayesian
  approach, negative binomial CUSUM methods and a detector based on
  generalized likelihood ratios.
  Furthermore, inference methods for the retrospective infectious disease
  model in Held et al (2005), Held et al (2006), Paul et al (2008) and
  Paul and Held (2011) are provided. A novel CUSUM approach combining
  logistic and multinomial logistic modelling is also included.
  Continuous self-exciting spatio-temporal point processes are modeled
  through additive-multiplicative conditional intensities as described in
  H\enc{ö}{oe}hle (2009) ("twinSIR", discrete space) and Meyer et al (2012)
  ("twinstim", continuous space).
  The package contains several real-world datasets, the ability to
  simulate outbreak data, visualize the results of the monitoring in
  temporal, spatial or spatio-temporal fashion.
} 

\details{
  \tabular{ll}{
    Package: \tab \Sexpr[stage=build]{meta$Package}\cr
    Version: \tab \Sexpr[stage=build]{meta$Version}\cr
    License: \tab \Sexpr[stage=build]{meta$License}\cr
    URL:     \tab \url{http://surveillance.r-forge.r-project.org/}\cr
    % cannot use \Sexpr within \url; and results=rd with "\\url"-paste is buggy
  }
  
  \pkg{surveillance} is an \R package implementing statistical methods
  for the retrospective modeling and prospective change-point detection
  in time series of counts, proportions and categorical data. The main
  application is in the detection of aberrations in routine collected
  public health data seen as univariate and multivariate time series of
  counts or point-processes. However, applications could just as well
  originate from environmetrics, econometrics or social sciences. As
  many methods rely on statistical process control methodology, the
  package is thus also relevant to quality control and reliability
  engineering.

  The fundamental data structure of the package is an S4 class
  \code{sts} wrapping observations, monitoring results and date handling
  for multivariate time series. Currently the package contains
  implementations typical outbreak detection procedures such as Stroup et
  al. (1989), Farrington et al., (1996), Rossi et al. (1999), Rogerson
  and Yamada (2001), a Bayesian approach (H\enc{ö}{oe}hle, 2007),
  negative binomial CUSUM methods (H\enc{ö}{oe}hle and Mazick, 2009), and a
  detector based on generalized likelihood ratios (H\enc{ö}{oe}hle
  and Paul, 2008). However, also CUSUMs for the prospective change-point
  detection in binomial, beta-binomial and multinomial time series is
  covered based on generalized linear modelling.  This includes
  e.g. paired binary CUSUM described by Steiner et al. (1999) or paired
  comparison Bradley-Terry modelling described in H\enc{ö}{oe}hle
  (2010). The package contains several real-world datasets, the ability
  to simulate outbreak data, visualize the results of the monitoring in
  temporal, spatial or spatio-temporal fashion.

  Furthermore, the package contains inference methods for the
  retrospective infectious disease model in Held et al. (2005), Paul et
  al. (2008) ("algo.hhh") and Paul and Held (2011) ("hhh4") handling
  multivariate time series of counts. Furthermore, the fully Bayesian
  approach for univariate time series of counts from Held et al. (2006)
  ("twins") is also implemented. Self-exciting spatio-temporal point
  processes are modeled through additive-multiplicative conditional
  intensities as described in H\enc{ö}{oe}hle (2009) ("twinSIR") and
  Meyer et al (2012) ("twinstim").

  Altogether, the package allows the modelling and monitoring of
  epidemic phenomena in temporal and spatio-temporal contexts.
}

%% Author information is dynamically extracted from the DESCRIPTION file
\author{
  \Sexpr[stage=build]{authors <- unname(eval(parse(text=meta$"Authors@R")))}
  \Sexpr[stage=build]{formatPerson <- function(person, sort=FALSE) paste0(format(if (sort && length(person) > 1) person[order(unlist(person$family))] else person, include=c("given", "family")), collapse=", ")}
  \Sexpr[stage=build]{formatPerson(authors[grep("aut", authors$role)])}

  \Sexpr[stage=build]{maintainer <- authors[grep("cre", authors$role)]}
  Maintainer: \Sexpr[stage=build]{formatPerson(maintainer)} \email{\Sexpr[stage=build]{maintainer$email}}
}

%% Dynamically extract contributors from the DESCRIPTION file
%% and persons from inst/THANKS for acknowledgement:
\section{Acknowledgements}{
  Substantial contributions of code by:
  \Sexpr[stage=build]{contributors <- authors[grepl("ctb", authors$role) & !sapply(authors$family, is.null)]}
  \Sexpr[stage=build]{formatPerson(contributors, sort=TRUE)}.

  Furthermore, the authors would like to thank the following people
  for ideas, discussions, testing and feedback:
  \Sexpr[stage=build]{THANKSfile <- file(system.file("THANKS", package="surveillance", mustWork=TRUE), encoding="latin1")} % this re-encodes from "latin1" to the current native encoding
  \Sexpr[stage=build]{formatPerson(as.person(grep("^(#|[[:blank:]]*$)", readLines(THANKSfile), invert=TRUE, value=TRUE)), sort=TRUE)}.
  \Sexpr[stage=build]{close(THANKSfile)}
}

\references{
  See \code{citation(package="surveillance")}.
}

\keyword{ package }

\examples{
#Code from an early survey article about the package: Hoehle (2007)
#available from http://surveillance.r-forge.r-project.org/
\dontrun{demo(cost)}

#Code from a more recent book chapter about using the package for the
#monitoring of Danish mortality data (Hoehle and Mazick, 2010).
\dontrun{demo(biosurvbook)}
}
