\name{m1}
\alias{m1}
\alias{h1_nrwrp}
\alias{k1}
\alias{m2}
\alias{m3}
\alias{m4}
\alias{m5}
\alias{n1}
\alias{n2}
\alias{q1_nrwh}
\alias{q2}
\alias{s1}
\alias{s2}
\alias{s3}
\docType{data}
\encoding{latin1}

\title{RKI SurvStat Data}
\description{
  14 datasets for different diseases beginning in
  2001 to the 3rd Quarter of 2004 including their defined outbreaks.
  \itemize{
    \item \code{m1} 'Masern' in the 'Landkreis Nordfriesland' (Germany, Schleswig-Holstein)
    \item \code{m2} 'Masern' in the 'Stadt- und Landkreis Coburg' (Germany, Bayern)
    \item \code{m3} 'Masern' in the 'Kreis Leer' (Germany, Niedersachsen)
    \item \code{m4} 'Masern' in the 'Stadt- und Landkreis Aachen' (Germany, Nordrhein-Westfalen)
    \item \code{m5} 'Masern' in the 'Stadt Verden' (Germany, Niedersachsen)
    \item \code{q1\_nrwh} 'Q-Fieber' in the 'Hochsauerlandkreis' (Germany, Westfalen)
    and in the 'Landkreis Waldeck-Frankenberg' (Germany, Hessen)
    \item \code{q2} 'Q-Fieber' in '\enc{München}{Muenchen}' (Germany, Bayern)
    \item \code{s1} 'Salmonella Oranienburg' in Germany
    \item \code{s2} 'Salmonella Agona' in 12 'Bundesl\enc{ä}{ae}ndern' of Germany
    \item \code{s3} 'Salmonella Anatum' in Germany
    \item \code{k1} 'Kryptosporidiose' in Germany, 'Baden-W\enc{ü}{ue}rttemberg'
    \item \code{n1} 'Norovirus' in 'Stadtkreis Berlin Mitte' (Germany, Berlin)
    \item \code{n2} 'Norovirus' in 'Torgau-Oschatz' (Germany, Sachsen)
    \item \code{h1\_nrwrp} 'Hepatitis A' in 'Oberbergischer Kreis, Olpe, Rhein-Sieg-kreis'
    (Germany, Nordrhein-Westfalen) and 'Siegenwittgenstein Altenkirchen' (Germany, Rheinland-Pfalz)
  }
}
\usage{data(m1)}
\format{
  \code{disProg} objects each containing 209 observations (weekly on 52 weeks)
  \describe{
    \item{observed}{Number of counts in the corresponding week}
    \item{state}{Boolean whether there was an outbreak.}
  }
}
\source{
  Robert Koch-Institut: SurvStat:
  \url{https://survstat.rki.de/};
  m1 and m3 were queried on 10 November 2004. The rest during September 2004.
}
\examples{
data(k1)
survResObj <- algo.rki1(k1, control=list(range=27:192))
plot(survResObj, "RKI 1", "k1", firstweek=27, startyear=2002)
}
\keyword{datasets}
