################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Simple wrapper around functionality of the numDeriv and maxLik packages
### to check the score vector and the Fisher information matrix
### CAVE: the return values of both wrappers are not unified
###
### Copyright (C) 2012 Sebastian Meyer
### $Revision: 463 $
### $Date: 2012-12-06 17:26:39 +0100 (Do, 06. Dez 2012) $
################################################################################


checkDerivatives.numDeriv <- function(ll, score, fisher, par,
                                      method="Richardson",
                                      method.args=list(), ...)
{
    cat("Checking analytical score vector using numDeriv::grad() ...\n")
    nsc <- numDeriv::grad(ll, par, method, method.args, ...)
    asc <- score(par, ...)
    print(all.equal(asc, nsc, check.attributes=FALSE))
    cat("Checking analytical Fisher information matrix using numDeriv::hessian() ...\n")
    if (length(par) > 50)
        cat("NOTE: this might take several minutes considering length(par) =",
            length(par), "\n")
    nfi <- -numDeriv::hessian(ll, par, "Richardson", method.args, ...)
    afi <- fisher(par, ...)
    print(all.equal(afi, nfi, check.attributes=FALSE))
    invisible(list(score  = list(analytic=asc, numeric=nsc),
                   fisher = list(analytic=afi, numeric=nfi)))
}


checkDerivatives.maxLik <- function(ll, score, fisher, par, eps=1e-6,
                                    print=FALSE, ...)
{
    cat("Checking analytical score and Fisher using maxLik::compareDerivatives() ...\n")
    res <- maxLik::compareDerivatives(
                   f=ll, grad=score,
                   hess=function (par, ...) -fisher(par, ...),
                   t0=par, eps=eps, print=print, ...)
    cat("Comparison of score vectors:\n")
    print(all.equal(res$compareGrad$analytic, drop(res$compareGrad$numeric),
                    check.attributes=FALSE))
    cat("Comparison of Fisher information matrices:\n")
    print(all.equal(res$compareHessian$analytic, drop(res$compareHessian$numeric),
                    check.attributes=FALSE))
    invisible(res)
}
