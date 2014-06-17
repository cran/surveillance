################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Permutation test to compare the means of paired samples
###
### Copyright (C) 2010-2012 Michaela Paul
### $Revision: 926 $
### $Date: 2014-05-22 13:47:29 +0200 (Thu, 22 May 2014) $
################################################################################

permutationTest <- function(score1, score2, nPermutation=9999, plot=FALSE)
{
    stopifnot((nTime <- length(score1)) == length(score2),
              !is.na(score1), !is.na(score2))
    
    meanScore1 <- mean(score1)
    meanScore2 <- mean(score2)
    diffObserved <- meanScore1 - meanScore2
    
    diffMean <- replicate(nPermutation, {
        sel <- rbinom(nTime, size=1, prob=0.5)
        g1 <- (sum(score1[sel==0]) + sum(score2[sel==1]))/nTime
        g2 <- (sum(score1[sel==1]) + sum(score2[sel==0]))/nTime
        g1 - g2
    })

    if (plot) {
        hist(diffMean, nclass=50, prob=TRUE,xlab="Difference between means",main="")
        abline(v=diffObserved,col=4)
    }
    
    pVal <- (1+sum(abs(diffMean)>=abs(diffObserved))) / (nPermutation+1)

    pTtest <- t.test(score1, score2, paired=TRUE)$p.value
    
    cat("mean difference=", diffObserved,
        "\tp(permutation) =", pVal,
        "\tp(paired t-test) =", pTtest, "\n")
    
    list(diffObs=diffObserved, pVal.permut=pVal, pVal.t=pTtest)
}
