/*******************************************************************************
// Determine potential triggering events close in space and time
//
// Copyright (C) 2016,2021,2024 Sebastian Meyer
// 
// This file is part of the R package "surveillance",
// free software under the terms of the GNU General Public License, version 2,
// a copy of which is available at https://www.R-project.org/Licenses/.
*******************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <math.h>

SEXP determineSources(
    SEXP eventTimesSEXP, SEXP eps_tSEXP,
    SEXP eventCoordsSEXP, SEXP eps_sSEXP,
    SEXP eventTypes0SEXP, SEXP qmatrixSEXP, SEXP nTypesSEXP
){
    int N = LENGTH(eventTimesSEXP), i;

    // Get pointers to input vectors
    double *eventTimes = REAL(eventTimesSEXP);
    double *eps_t = REAL(eps_tSEXP);
    double *xcoords = REAL(eventCoordsSEXP);
    double *ycoords = xcoords + N; // move the pointer to the second column
    double *eps_s = REAL(eps_sSEXP);
    int *eventTypes0 = INTEGER(eventTypes0SEXP); // 0-based type index
    int *qmatrix = LOGICAL(qmatrixSEXP);
    int nTypes = asInteger(nTypesSEXP);

    // Allocate output list of sources for each event
    SEXP sourcesSEXP = PROTECT(allocVector(VECSXP, N));
    
    for (i = 0; i < N; i++) {
        SEXP sourceVec = PROTECT(allocVector(INTSXP, N)); // maximum length
        int nSources = 0;
        // Check for each event j if it could be a source of i
        for (int j = 0; j < N; j++) {
            if (// 1. type of j can infect type of i
                qmatrix[eventTypes0[j] + eventTypes0[i] * nTypes] &&
                // 2. j is still infectious when i happens
                // (use "<" not "<=": CIF is left-continuous, no self-infection)
                (eventTimes[j] < eventTimes[i]) &&
                (eventTimes[j] + eps_t[j] >= eventTimes[i]) &&
                // 3. j reaches i
                // (use "hypot" for Euclidean distance between j and i)
                hypot(xcoords[j] - xcoords[i], ycoords[j] - ycoords[i]) <= eps_s[j])
            {
                INTEGER(sourceVec)[nSources] = j + 1; // 1-based index for R
                nSources++;
            }
        }
        sourceVec = lengthgets(sourceVec, nSources); // trim to actual length
        SET_VECTOR_ELT(sourcesSEXP, i, sourceVec);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return sourcesSEXP;
}
