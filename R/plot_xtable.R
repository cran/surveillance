###################################################
### chunk number 1: 
###################################################

compMatrix.writeTable <- function(compMatrix){
        require(xtable)
        compMatrix.table <- xtable(compMatrix, display = c("s", rep("d", 4), rep("G", 4)), vsep = "|")

        return(compMatrix.table)
}


