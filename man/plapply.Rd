\name{plapply}
\alias{plapply}
\title{Verbose and Parallel \code{lapply}}

\description{
Verbose and parallelized version of \code{lapply} wrapping around
\code{\link[parallel]{mclapply}} and \code{\link[parallel]{parLapply}}
in the base package \pkg{parallel}.
This wrapper can take care of the \code{.Random.seed} and
print progress information (not for cluster-based parallelization).
With the default arguments it equals \code{lapply}
enriched by a progress bar.
}

\usage{
plapply(X, FUN, ...,
        .parallel = 1, .seed = NULL, .verbose = TRUE)
}

\arguments{
  \item{X,FUN,\dots}{see \code{\link{lapply}}.}
  \item{.parallel}{
    the number of processes to use in parallel operation, or a
    \code{"cluster"} object (see \code{\link[parallel]{makeCluster}}).
    If a number, \code{\link[parallel]{mclapply}} (forking) is used on
    Unix-alikes, whereas on Windows \code{\link[parallel]{parLapply}} is
    used on a newly created cluster of the specified size, which is
    stopped when exiting the function.
    By default (\code{.parallel = 1}), the basic \code{\link{lapply}} is used.
  }
  \item{.seed}{
    If set (non-\code{NULL}), results involving random number generation
    become reproducible. If using a cluster (see the \code{.parallel}
    argument), \code{\link[parallel]{clusterSetRNGStream}} is called
    with the specified \code{.seed} before running \code{parLapply}.
    Otherwise, \code{\link{set.seed}(.seed)} is called and the
    \code{\link{RNGkind}} is changed to \code{"L'Ecuyer-CMRG"} if
    \code{.parallel > 1} (see the section on random numbers in the
    documentation of \code{mcparallel} in package \pkg{parallel}).
    % no link to mcparallel since it is not available on Windows (R-3.1.2)
    If \code{.seed} is non-\code{NULL},
    the original \code{\link{.Random.seed}} will be restored
    \code{on.exit} of the function.
  }
  \item{.verbose}{
    if and how progress information should be displayed, i.e., what to
    do on each exit of \code{FUN}. This is unsupported and ignored for
    cluster-based parallelization and primitive \code{FUN}ctions.
    The default (\code{TRUE}) will show a \code{\link{txtProgressBar}}
    (if \code{.parallel = 1} in an \code{\link{interactive}} \R session)
    or \code{cat(".")} (otherwise). Other choices for the dot are possible
    by specifying the desired symbol directly as the \code{.verbose}
    argument. Alternatively, \code{.verbose} may be any custom call or
    expression to be executed \code{\link{on.exit}} of \code{FUN} and
    may thus involve any objects from the local evaluation environment.
  }
}

\value{
  a list of the results of calling \code{FUN} on each value of \code{X}.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link[parallel]{mclapply}} and \code{\link[parallel]{parLapply}}
}

\examples{
## example inspired by help("lapply")
x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))

## if neither parallel nor verbose then this simply equals lapply()
plapply(x, quantile, probs = 1:3/4, .verbose = FALSE)

## verbose lapply() -- not really useful for such fast computations
res <- plapply(x, quantile, probs = 1:3/4, .verbose = TRUE)
res <- plapply(x, quantile, probs = 1:3/4, .verbose = "|")
res <- plapply(x, quantile, probs = 1:3/4,
               .verbose = quote(cat("length(x) =", length(x), "\n")))

## setting the seed for reproducibility of results involving the RNG
samp <- plapply(as.list(1:3), runif, .seed = 1)

## parallel lapply()
res <- plapply(x, quantile, probs = 1:3/4, .parallel = 2)

## using a predefined cluster
library("parallel")
cl <- makeCluster(getOption("cl.cores", 2))
res <- plapply(x, quantile, probs = 1:3/4, .parallel = cl)
stopCluster(cl)
}

\keyword{iteration}
\keyword{list}
