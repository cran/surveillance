################################################################################
### Yet another P-value formatter, using R's format.pval()
###
### Copyright (C) 2013,2015,2017 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

formatPval <- function (pv, eps = 1e-4, scientific = FALSE, ...)
{
    format1 <- function (p)
        format.pval(p, digits = if (p < 10*eps) 1 else 2, eps = eps,
                    nsmall = 2, scientific = scientific, ...)
    vapply(X = pv, FUN = format1, FUN.VALUE = "", USE.NAMES = TRUE)
}
