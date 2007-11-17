# -------------  class sts  ----------------------------------------

setClass( "sts", representation(week = "numeric",
                                freq = "numeric",
                                start = "numeric",
                                observed = "matrix",
                                state = "matrix",
                                alarm = "matrix",
                                upperbound  = "matrix",
                                neighbourhood= "matrix",
                                populationFrac= "matrix",
#                                lvl = "vector",
                                map = "SpatialPolygonsDataFrame",
                                control = "list"))

