\name{bestCombination}
\alias{bestCombination}
\encoding{latin1}

\title{Partition of a number into two factors}
\description{
       Given a prime number factorization \code{x}, \code{best combination}
       partitions \code{x} into two groups, such that the product of the numbers
       in group one is as similar as possible to the product
       of the numbers of group two. This is useful in magic.dim 
}
\usage{
      bestCombination(x)
}
\arguments{
      \item{x}{prime number factorization}
}
\value{Returns a vector c(prod(set1),prod(set2))
      }
%\keyword{}
