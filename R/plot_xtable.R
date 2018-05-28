###################################################
### chunk number 1:
###################################################

compMatrix.writeTable <- function (compMatrix)
{
    .Deprecated(package = "surveillance")
    xtable(compMatrix, display = c("s", rep("d", 4), rep("G", 4)), vsep = "|")
}
