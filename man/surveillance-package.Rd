% Process \Sexpr code at build time; like bib macros (not yet supported by pkgdown).
\RdOpts{stage=build}
% NOTE: in R 2.12.0 - 4.0.5, \RdOpts{} had no effect (PR#18073)

\encoding{UTF-8}
\name{surveillance-package}
\alias{surveillance-package}
\alias{surveillance}
\docType{package}

\title{\pkg{surveillance}: \packageTitle{surveillance}}

\description{
  The \R package \pkg{surveillance} implements statistical methods for the
  retrospective modeling and prospective monitoring of epidemic phenomena
  in temporal and spatio-temporal contexts.
  Focus is on (routinely collected) public health surveillance data,
  but the methods just as well apply to data from environmetrics,
  econometrics or the social sciences. As many of the monitoring methods
  rely on statistical process control methodology, the package is
  also relevant to quality control and reliability engineering.
} 

\details{
  The package implements many typical outbreak detection procedures such
  as Stroup et al. (1989), Farrington et al. (1996), Rossi et al. (1999),
  Rogerson and Yamada (2001), a Bayesian approach (Höhle, 2007),
  negative binomial CUSUM methods (Höhle and Mazick, 2009), and a
  detector based on generalized likelihood ratios (Höhle and Paul, 2008),
  see \code{\link{wrap.algo}}.
  Also CUSUMs for the prospective change-point detection in binomial,
  beta-binomial and multinomial time series are covered based on
  generalized linear modeling, see \code{\link{categoricalCUSUM}}.
  This includes, e.g., paired comparison Bradley-Terry modeling described
  in Höhle (2010), or paired binary CUSUM
  (\code{\link{pairedbinCUSUM}}) described by Steiner et al. (1999).
  The package contains several real-world datasets, the ability
  to simulate outbreak data, visualize the results of the monitoring in
  temporal, spatial or spatio-temporal fashion. In dealing with time
  series data, the fundamental data structure of the package is the S4
  class \code{\link{sts}} wrapping observations, monitoring results and
  date handling for multivariate time series.
  A recent overview of the available monitoring procedures is
  given by Salmon et al. (2016).

  For the retrospective analysis of epidemic spread, the package
  provides three endemic-epidemic modeling frameworks with
  tools for visualization, likelihood inference, and simulation.
  The function \code{\link{hhh4}} offers inference methods for the
  (multivariate) count time series model of Held et al. (2005),
  including some of its extensions (Paul and Held, 2011; Meyer and
  Held, 2014). See \vignette{hhh4} for a general introduction
  and \vignette{hhh4_spacetime} for a discussion and
  illustration of spatial \code{hhh4} models.
  Self-exciting point processes are modeled through endemic-epidemic
  conditional intensity functions.
  \code{\link{twinSIR}} (Höhle, 2009) models the
  susceptible-infectious-recovered (SIR) event history of a 
  fixed population, e.g, epidemics across farms or networks;
  see \vignette{twinSIR} for an illustration.
  \code{\link{twinstim}} (Meyer et al., 2012) fits spatio-temporal point
  process models to point patterns of infective events, e.g.,
  time-stamped geo-referenced surveillance data on infectious disease
  occurrence; see \vignette{twinstim} for an illustration.
  A recent overview of the implemented space-time modeling frameworks
  for epidemic phenomena is given by Meyer et al. (2017).
}

%% Author information is dynamically extracted from the DESCRIPTION file
\author{
  % \packageAuthor{surveillance} would simply give the Author field text,
  % including all roles and without markup, so:
  \Sexpr[results=rd]{
    pkgdir <- tools:::Rd_macros_package_dir()  # support R CMD Rd2pdf pkg
    aartxt <- read.dcf(file.path(pkgdir, "DESCRIPTION"), "Authors@R")
    aar <- eval(str2lang(aartxt))
    authors <- aar[grep("aut", aar$role)]
    # FIXME: why is there is no toRd.person() or style="Rd" for format.person() ?
    authors$comment <- lapply(authors$comment, function (x) {
      if (hasName(x, "ORCID"))
        sprintf("ORCID: \\\\href{https://orcid.org/\%s}{\%s}", x["ORCID"], x["ORCID"])
      else ""
    })
    c("\\\\itemize{", paste("\\\\item",
        format(authors, include = c("given", "family", "comment"))
      ), "}")
  }

  % \packageMaintainer{surveillance} lacks \email markup, just like
  %   format(aar[grep("cre", aar$role)], include=c("given","family","email"))
  % which does not support style="Rd", so:
  Maintainer: Sebastian Meyer \email{seb.meyer@fau.de}

  %% Dynamically extract contributors from the DESCRIPTION file
  %% and persons from inst/THANKS for acknowledgement:
  The authors acknowledge code contributions by
  \Sexpr{
    contributors <- aar[grepl("ctb", aar$role) & !sapply(aar$family, is.null)]
    paste0(format(contributors, include = c("given", "family")), collapse = ", ")
  }.

  Furthermore, the authors would like to thank the following people
  for ideas, discussions, testing and feedback:
  \Sexpr{
    find_inst_file <- function (file) { # support R CMD Rd2pdf in source package
      if(dir.exists(file.path(pkgdir, "inst")))
        file.path(pkgdir, "inst", file) else
        file.path(pkgdir, file)
    }
    ## internally, Rd objects are assumed to be in UTF-8, see parse_Rd()
    thanks <- readLines(find_inst_file("THANKS"), encoding = "UTF-8")
    toString(grep("^(#|$)", trimws(thanks), invert=TRUE, value=TRUE))
  }.
}

\section{\file{CITATION}}{
  \Sexpr[results=hide]{
    CITATION <- readCitationFile(find_inst_file("CITATION"))  # gives marked UTF-8 strings
  }
  To cite the \R package \pkg{surveillance} in publications, please use:
  \describe{
    \item{For the monitoring functions:}{%\bibshow{salmon.etal2016}
      \Sexpr[results=rd]{tools::toRd(CITATION[grep("Salmon", CITATION$author)])}
    }
    \item{For the modeling frameworks:}{%\bibshow{meyer.etal2017}
      \Sexpr[results=rd]{tools::toRd(CITATION[grep("Meyer", CITATION$author)])}
    }
  }
  See also \code{citation("surveillance")}.
}

\references{
  References to the underlying methodological papers are given on
  the specific help pages; some are also available in an
  installed BibTeX file:
\preformatted{
    system.file("REFERENCES.bib", package = "surveillance")
}
  %\bibshow{*}
}

\seealso{
  \url{https://surveillance.R-forge.R-project.org/}
}

\keyword{ package }

\examples{
## Additional documentation and illustrations of the methods are
## available in the form of package vignettes and demo scripts:
vignette(package = "surveillance")
demo(package = "surveillance")
}
