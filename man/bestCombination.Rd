\name{bestCombination}
\alias{bestCombination}
\title{Partition of a number into two factors}
\description{
  Given a prime number factorization \code{x}, \code{bestCombination}
  partitions \code{x} into two groups, such that the product of the numbers
  in group one is as similar as possible to the product
  of the numbers of group two. This is useful in \code{\link{magic.dim}}.
}
\usage{
bestCombination(x)
}
\arguments{
  \item{x}{prime number factorization}
}
\value{a vector \code{c(prod(set1),prod(set2))}}
\keyword{dplot}
