/*******************************************************************************
// Determine potential triggering events close in space and time
//
// Copyright (C) 2016,2021 Sebastian Meyer
// 
// This file is part of the R package "surveillance",
// free software under the terms of the GNU General Public License, version 2,
// a copy of which is available at https://www.R-project.org/Licenses/.
*******************************************************************************/

#include <Rcpp.h>

using namespace Rcpp;


// Euclidean distance of a set of points to a single point (x0, y0)
NumericVector distsN1(NumericVector x, NumericVector y, double x0, double y0)
{
        // hypot(x, y) is not (yet) vectorized by Rcpp sugar
        return sqrt(pow(x - x0, 2.0) + pow(y - y0, 2.0));
}


RcppExport SEXP determineSources(
    SEXP eventTimesSEXP, SEXP eps_tSEXP,
    SEXP eventCoordsSEXP, SEXP eps_sSEXP,
    SEXP eventTypesSEXP, SEXP qmatrixSEXP
){
    BEGIN_RCPP

    NumericVector eventTimes(eventTimesSEXP);
    NumericVector eps_t(eps_tSEXP);
    NumericMatrix eventCoords(eventCoordsSEXP);
    NumericVector eps_s(eps_sSEXP);
    IntegerVector eventTypes(eventTypesSEXP);
    LogicalMatrix qmatrix(qmatrixSEXP);

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
        
        return wrap(sources);

    END_RCPP
}
