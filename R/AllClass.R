# -------------  class sts  ----------------------------------------

sts <- setClass( "sts", representation(epoch = "numeric",  ##this slot used to be named week
                                freq = "numeric",
                                start = "numeric",
                                observed = "matrix",
                                state = "matrix",
                                alarm = "matrix",
                                upperbound  = "matrix",
                                neighbourhood= "matrix",
                                populationFrac= "matrix",
                                map = "SpatialPolygons",
                                control = "list",
#New slots added to handle proportion time series 
                                epochAsDate="logical",
                                multinomialTS="logical"))



######################################################################
# Definition of the stsBP class for backprojections.
######################################################################

setClass( "stsBP", slots = list(ci = "array", lambda = "array"),
         contains = "sts")

######################################################################
# Definition of the stsNC class for nowcasts.
######################################################################

setClass( "stsNC", slots=list(reportingTriangle = "matrix", predPMF= "list", pi = "array", truth= "sts", delayCDF = "list", SR = "array"),
         contains = "sts")
