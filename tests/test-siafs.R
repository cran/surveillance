library("surveillance")

### check siaf derivatives and Fcircle
testsiaf <- function (siaf, pargrid, type=1)
{
    ints <- surveillance:::checksiaf.Fcircle(
        siaf$Fcircle, siaf$f, pargrid, type=type, rs=c(1,5,10,50,100), nGQ=15)
    stopifnot(isTRUE(all.equal(ints[,1], ints[,2], tolerance=0.005)))

    maxRelDiffs_deriv <- surveillance:::checksiaf.deriv(
        siaf$deriv, siaf$f, pargrid, type=type)
    stopifnot(maxRelDiffs_deriv < 0.0001)
}

### Gaussian kernel
testsiaf(siaf.gaussian(1, F.adaptive=TRUE), as.matrix(c(1,3,6)))

### Power law kernels
testsiaf(siaf.powerlaw(), t(c(-1.5, 0.5)))
testsiaf(siaf.powerlawL(), t(c(-1, 0.4)))
