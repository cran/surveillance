################################################################################
### Generate a color palette via the colorspace package
###
### Copyright (C) 2007 Michael Hoehle, 2012-2014,2017 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

hcl.colors <- function (ncolors=100, use.color=TRUE)
{
    GYR <- if (requireNamespace("colorspace", quietly=TRUE)) {
        ## the Zeil-ice colors
        colorspace::heat_hcl(ncolors, h=c(0,120),
                             c=if (use.color) c(90,30) else c(0,0),
                             l=c(50,90), power=c(0.75, 1.2))
    } else {
        if (use.color) heat.colors(ncolors) else grey.colors(ncolors)
    }
    return(rev(GYR))
}
