################################################################################
### Conditional lapply
###
### Copyright (C) 2012,2017 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

### clapply uses lapply if X is a list and otherwise applies FUN directly to X.
### The result is always a list (of length 1 in the latter case).

clapply <- function (X, FUN, ...)
{
    if (is.list(X)) lapply(X, FUN, ...) else list(FUN(X, ...))
}
