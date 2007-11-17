\name{find.kh}
\alias{find.kh}
\title{Determine the k and h values in a standard normal setting}
\description{
        Given a specification of the average run length in the (a)cceptance
        and (r)ejected setting determine the k and h values in a standard
        normal setting.
}
\usage{
  find.kh(ARLa = 500, ARLr = 7, sided = "one", verbose=FALSE)
}
\arguments{
        \item{ARLa}{average run length in acceptance setting, aka. in control state. Specifies the number of observations before false alarm.}
        \item{ARLr}{average run length in rejection state, aka. out of control state. Specifies the number of observations before an increase is detected (i.e. detection delay)}
        \item{sided}{one-sided cusum scheme}
        \item{verbose}{gives extra information about the root finding process}
}
\value{
     Returns a list with reference value k and decision interval h.
}
\details{
        Functions from the spc package are used in a simple univariate
        root finding problem.
}
%\seealso{}
%\author{}
\examples{

    find.kh(ARLa=500,ARLr=7,sided="one")
    find.kh(ARLa=500,ARLr=3,sided="one")
}
%\keyword{}
%\source{}

