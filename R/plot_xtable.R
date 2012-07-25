###################################################
### chunk number 1: 
###################################################

compMatrix.writeTable <- function (compMatrix)
{
    xtable::xtable(compMatrix, display = c("s", rep("d", 4), rep("G", 4)), vsep = "|")
}


