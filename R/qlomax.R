################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Simple implementation of the quantile function of the Lomax distribution
### (we could also use VGAM::qlomax, but this would be slightly slower)
###
### Copyright (C) 2012-2013 Sebastian Meyer
### $Revision: 2124 $
### $Date: 2018-05-11 10:10:31 +0200 (Fri, 11. May 2018) $
################################################################################

qlomax <- function (p, scale, shape) {
    .Deprecated("VGAM::qlomax", package = "surveillance")
    scale * ((1-p)^(-1/shape) - 1)
}
