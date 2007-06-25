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
                                map = "SpatialPolygonsDataFrame",
                                control = "list"))

