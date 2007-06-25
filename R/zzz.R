###################################################
### chunk number 1: 
###################################################

.First.lib <- function(libname, pkgname) {
  #Load the necessary packages
  library(spc)
  library(maptools)

  #Load the CIdata thing
  data(CIdata, package=pkgname)

  #Read the table of the hypgeom_2F1 function for parameters c(1/3,2/3) and
  #5/3 -- atm this is computed for the values seq(0,10,by=0.01) and 11:100
  #Load the pre-evaluated Hypergeometric function for computing Anscombe residuals
  surveillance.gvar.hyp <- scan(file.path(.path.package('surveillance'),'data',"hypGeomSmall.txt"),quiet=TRUE)
  surveillance.gvar.z <- - c(0:1000/100, 11:100)

  #Load the C code library for the glr stuff
  library.dynam("surveillance", pkgname, libname)
}




