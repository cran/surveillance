/*******************************************************************************
// Determine potential triggering events close in space and time.
// Copyright (C) 2016  Sebastian Meyer <sebastian.meyer@ifspm.uzh.ch>
// 
// This program is part of the surveillance package,
// http://surveillance.r-forge.r-project.org,
// free software under the terms of the GNU General Public License, version 2,
// a copy of which is available at http://www.r-project.org/Licenses/.
*******************************************************************************/


#include <Rcpp.h>
using namespace Rcpp;


// Euclidean distance of a set of points to a single point (x0, y0)
NumericVector distsN1(NumericVector x, NumericVector y, double x0, double y0)
{
        // hypot(x, y) is not (yet) vectorized by Rcpp sugar
        return sqrt(pow(x - x0, 2.0) + pow(y - y0, 2.0));
}


// [[Rcpp::export]]
List determineSourcesC(
        NumericVector eventTimes, NumericVector eps_t,
        NumericMatrix eventCoords, NumericVector eps_s,
        IntegerVector eventTypes, LogicalMatrix qmatrix
){
        int N = eventTimes.size();
        NumericVector removalTimes = eventTimes + eps_t;
        NumericMatrix::Column xcoords = eventCoords(_,0);
        NumericMatrix::Column ycoords = eventCoords(_,1);
        List sources(N);
        LogicalVector infectivity(N);
        LogicalVector proximity(N);
        LogicalVector matchType(N);
        LogicalVector typeInfective(qmatrix.nrow());
        IntegerVector eventTypes0 = eventTypes - 1; // for correct indexing
        IntegerVector idx = seq_len(N);
        
        for (int i = 0; i < N; ++i) {
                infectivity = (eventTimes < eventTimes[i]) &
                        (removalTimes >= eventTimes[i]);
                // "<" not "<=" because CIF is left-continuous.
                // Also guarantees no self-infection.
                proximity = distsN1(xcoords, ycoords, eventCoords(i,0), eventCoords(i,1)) <= eps_s;
                typeInfective = qmatrix(_,eventTypes0[i]);
                //<- logical vector indicating for each type if it could infect type of i
                matchType = typeInfective[eventTypes0];

                sources[i] = idx[infectivity & proximity & matchType];
        }
        
        return sources;
}



// The following R code will be run automatically after compilation by
// Rcpp::sourceCpp("~/Projekte/surveillance/pkg/src/determineSources.cc")

/*** R
data("imdepi", package="surveillance")
sources <- imdepi$events$.sources
tail(sources)

eventTimes <- imdepi$events$time
eps.t <- imdepi$events$eps.t
eventCoords <- coordinates(imdepi$events)
eps.s <- imdepi$events$eps.s
eventTypes <- imdepi$events$type
qmatrix <- imdepi$qmatrix

sourcesC <- determineSourcesC(eventTimes, eps.t, eventCoords, eps.s, as.integer(eventTypes), qmatrix)
tail(sourcesC)
stopifnot(identical(sources, sourcesC))

library("microbenchmark")
microbenchmark(
    determineSourcesC(eventTimes, eps.t, eventCoords, eps.s, as.integer(eventTypes), qmatrix),
    surveillance:::determineSources.epidataCS(imdepi, method = "R"),
    times = 50)
*/



/*** This is how tedious the function would look like without Rcpp attributes:
RcppExport SEXP determineSourcesCSEXP(SEXP eventTimesSEXP, SEXP eps_tSEXP,
                                      SEXP eventCoordsSEXP, SEXP eps_sSEXP,
                                      SEXP eventTypesSEXP, SEXP qmatrixSEXP)
{
        NumericVector eventTimes(eventTimesSEXP);
        NumericVector eps_t(eps_tSEXP);
        NumericMatrix eventCoords(eventCoordsSEXP);
        NumericVector eps_s(eps_sSEXP);
        IntegerVector eventTypes(eventTypesSEXP);
        LogicalMatrix qmatrix(qmatrixSEXP);
        
[... insert body of the above determineSourcesC here but replace return statement by ...]

        return wrap(sources);
}
*/
