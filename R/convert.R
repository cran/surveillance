######################################################################
# Convert old STS class with slotname "week" to new sts
# where this slot is called "epoch"
#
# Parameters:
#  oldSTS - an "old school" sts S4 object (having a slot called "week")
#
# Returns:
#  a "new school" STS object having a slot called "epoch"
######################################################################

convertSTS <- function(oldSTS) {#.Object, epoch, start=c(2000,1), freq=52, observed, state=0*observed, map=NULL, neighbourhood=NULL, populationFrac=NULL,alarm=NULL,upperbound=NULL, control=NULL,epochAsDate=FALSE,multinomialTS=FALSE) {
  new("sts",
      epoch=oldSTS@week, #here is the action,
      freq=oldSTS@freq,
      start=oldSTS@start,
      observed=oldSTS@observed,
      state=oldSTS@state,
      alarm=oldSTS@alarm,
      upperbound=oldSTS@upperbound,
      neighbourhood=oldSTS@neighbourhood,
      populationFrac=oldSTS@populationFrac,
      map=oldSTS@map,
      control=oldSTS@control,
      epochAsDate=oldSTS@epochAsDate,
      multinomialTS=oldSTS@multinomialTS)
}

doIt <- function() {
  #Load library
  library("surveillance")


  ######################################################################
  #Convert all sts S4 objects in the data directory
  ######################################################################
  
  #Momo
  data("momo")
  momo <- convertSTS(momo)
  save(file="data/momo.RData",list=c("momo"))

  #Flu data from Michaela
  load(file="data/flu-BYBW.RData")
  sts.flu <- convertSTS(sts.flu)
  save(file="data/flu-BYBW.RData",list=c("map","sts.flu","wji"))
  
  #Measles data
  data("measlesDE")
  measles.land <- convertSTS(measles.land)
  save(file="measlesDE.RData",list=c("measles.land","vac2006"))

  #Abattoir
  data("abattoir")
  abattoir <- convertSTS(abattoir)
  save(file="abattoir.RData",list=c("abattoir"))

  #Deleval
  data("deleval")
  deleval <- convertSTS(deleval)
  save(file="deleval.RData",list=c("deleval"))
  
}
