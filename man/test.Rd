\name{test}
\alias{test}
\encoding{latin1}
\title{Print xtable for several diseases and the summary}
\description{Just a test method}
\usage{
test(data = c("k1", "m5"), range = 157:339)
}
\arguments{
  \item{data}{vector of abbreviations for the diseases}
  \item{range}{timepoints to evaluate}
}
\value{
  printed LaTeX xtables
}
\details{
  The specified datasets are read, corrected, enlarged and sent to
  the RKI 1, RKI 2, RKI 3 and Bayes system. The quality values are computed
  and printed for each disease as latex table. Additionally a summary latex table
  for all diseases is printed.
}
\author{M. \enc{Höhle}{Hoehle}, A. Riebler, C. Lang}
\examples{
test(
    c("m1", "m2", "m3", "m4", "m5", "q1_nrwh", "q2",
      "s1", "s2", "s3", "k1", "n1", "n2", "h1_nrwrp")
    )
}
\keyword{misc}
