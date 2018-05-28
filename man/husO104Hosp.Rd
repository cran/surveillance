\encoding{latin1}
\name{husO104Hosp}
\alias{husO104Hosp}
\docType{data}
\title{Hospitalization date for HUS cases of the STEC outbreak in Germany, 2011}
\description{
  Data contain the date of hospitalization for 630 hemolytic-uremic
  syndrome (HUS) cases during the large STEC outbreak in Germany,
  2011. Note: Only HUS cases which ultimately had a hospitalization
  date available/reported are included in the data set. The total number
  of HUS cases during the outbreak was 855 -- see \enc{Höhle}{Hoehle} and
  an der Heiden (2014) as well as Frank et al. (2011) for details.
  
  For each HUS case the attribute \code{dHosp} contains the date of
  hospitalization and the attribute \code{dReport} contains the date of
  first arrival of this hospitalization date at the Robert Koch
  Institute (RKI). As described in \enc{Höhle}{Hoehle} and an der Heiden
  (2014) the mechanisms of the delay were complicated and should be
  interpreted with care. For example, the case report could have arrived
  earlier, but without information about the hospitalization date.

  The resulting reporting triangle corresponds to Fig. 1 of the Web
  appendix of \enc{Höhle}{Hoehle} and an der Heiden (2014). This means
  that the reports which arrived with a delay longer than 15 days are
  set to have have arrived after 15 days. Altogether, this gives small
  discrepancies when compared with the results of the paper. However, as
  mentioned in the paper, longer delays were not very relevant for the
  nowcasting.
}

\usage{data(husO104Hosp)}
\format{
  A \code{data.frame} object. 
}

\source{
  Data were collected during the outbreak as part of the mandatory
  reporting of notifiable diseases in Germany (Faensen et al.,
  2006). Here, reports are transmitted from the local health authorities
  via the state health authorities to the Robert Koch Institute, Berlin.
  The resulting reporting triangle corresponds to Fig. 1 of the Web
  appendix of \enc{Höhle}{Hoehle} and an der Heiden (2014).
}

\references{
  \enc{Höhle}{Hoehle} M and an der Heiden, M (2014). Bayesian Nowcasting
  during the STEC O104:H4 Outbreak in Germany, 2011, In revision for
  Biometrics.

  Frank C, Werber D, Cramer JP, Askar M, Faber M, an der Heiden M,
  Bernard H, Fruth A, Prager R, Spode A, Wadl M, Zoufaly A, Jordan S,
  Kemper MJ, Follin P, \enc{Müller}{Mueller} L, King LA, Rosner B,
  Buchholz U, Stark K, Krause G; HUS Investigation Team (2011).
  Epidemic Profile of Shiga-Toxin Producing Escherichia coli O104:H4
  Outbreak in Germany, N Engl J Med. 2011 Nov 10;365(19):1771-80.

  Faensen D, Claus H, Benzler J, Ammon A, Pfoch T, Breuer T, Krause G (2014).
  SurvNet@RKI - a multistate electronic reporting system for communicable
  diseases, Euro Surveillance, 2006;11(4):100-103.
}

\keyword{datasets}
