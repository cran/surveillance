################################################################################
### Check if an R object is scalar, i.e., a numeric vector of length 1
###
### Copyright (C) 2009,2017 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

isScalar <- function (x) {
    length(x) == 1L && is.vector(x, mode = "numeric")
}
