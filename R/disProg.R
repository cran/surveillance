###################################################
### chunk number 1:
###################################################
create.disProg <- function(week, observed, state, start=c(2001,1), freq=52, neighbourhood=NULL, populationFrac=NULL,epochAsDate=FALSE){
  ## issue a deprecation warning if not internally called
  if (!isTRUE(packageName(parent.frame()) == .packageName))
    .Deprecated("sts")
  
  namesObs <-colnames(observed)

  # check whether observed contains only numbers
  stopifnot(is.numeric(observed))

  #univariate timeseries ?
  if(is.vector(observed)){
    observed <- matrix(observed,ncol=1)
    namesObs <- "observed"
  } else {  # ensure we have a matrix
    observed <- as.matrix(observed)
  }

  if(missing(state)){
    state <- 0*observed
  } else if(is.vector(state)){
    state <- matrix(state,ncol=1)
  } else {
    state <- as.matrix(state)
  }

  #check number of columns of observed and state
  nAreas <- ncol(observed)
  nObs <- nrow(observed)
  if(ncol(state) != nAreas){
    #if there is only one state-vector for more than one area, repeat it
    if(ncol(state)==1) {
      state <- matrix(rep(state,nAreas),ncol=nAreas,byrow=FALSE)
    } else {
      stop("wrong dimensions of 'observed' and 'state'")
    }
  }

  #check neighbourhood matrix
  # neighbourhood can be a matrix or an array of dimension c(nAreas,nAreas, nrow(observed))
  if(!is.null(neighbourhood) ) {
    dimNhood <- dim(neighbourhood)
    if(!(length(dimNhood) %in% 2:3) ||
       any(dimNhood[1:2] != nAreas) ||
       (length(dimNhood)==3 && dimNhood[3] != nrow(observed))) {
      stop('wrong dimensions of neighbourhood matrix')
    }
  } else {
     # no neighbourhood specified
     neighbourhood <- matrix(NA,nrow=nAreas,ncol=nAreas)
  }

  if(is.null(populationFrac)) {
    populationFrac <- matrix(1/nAreas,nrow=nObs,ncol=nAreas)
  } else {
    # make sure populationFrac is a matrix
    populationFrac <- as.matrix(populationFrac)
    # check dimensions
    if(nrow(populationFrac)!= nObs | ncol(populationFrac)!= nAreas)
      stop("dimensions of 'populationFrac' and 'observed' do not match")
    # check whether populationFrac contains only numbers
    if(!is.numeric(populationFrac))
      stop("'populationFrac' must be a numeric matrix")
  }

  #labels for observed and state
  if(is.null(namesObs)){
    namesObs <- paste0("observed",1:nAreas)
  }

  colnames(observed) <- namesObs
  colnames(state) <- namesObs

  res <- list("week"=week, "observed"=observed, "state"=state, "start"=start, "freq"=freq,  "neighbourhood"=neighbourhood, "populationFrac"=populationFrac,"epochAsDate"=epochAsDate)
  class(res) <- "disProg"
  return(res)
}

print.disProg <- function(x, ...) {
  cat( "-- An object of class disProg -- \n" )
  cat( "freq:\t\t", x$freq,"\n" )
  cat( "start:\t\t", x$start,"\n" )
  cat( "dim(observed):\t", dim(x$observed), "\n\n")

  n <- 1
  cat("Head of observed:\n")
  print(head(x$observed,n))

  #cat("\nhead of neighbourhood:\n")
  #print( head(x$neighbourhood,n))
}


###################################################
### chunk number 3:
###################################################
aggregate.disProg <- function(x,...){
  #aggregate observed counts
  observed <- apply(x$observed,MARGIN=1,sum)
  #aggregate states
  state <- apply(x$state,MARGIN=1,sum)
  state[state > 1] <- 1

  #create univariate disProg object
  x <- create.disProg(week=x$week, observed=observed, state=state, freq=x$freq,start=x$start)
  return(x)
}


###################################################
### chunk number 4:
###################################################

## legacy code removed in surveillance 1.20.0
## now plotting via syntax transformation to stsplot_time(disProg2sts(x))
plot.disProg <- function(x, title = "", xaxis.years = TRUE,
                         startyear = x$start[1], firstweek = x$start[2],
                         as.one = TRUE, same.scale = TRUE, ...)
{
  cl <- match.call()
  cl[[1]] <- quote(surveillance::stsplot_time)
  stopifnot(!missing(x))
  cl$x <- substitute(surveillance::disProg2sts(x))
  names(cl)[names(cl) == "title"] <- "main"
  if (!xaxis.years) cl["xaxis.labelFormat"] <- list(NULL)
  cl$xaxis.years <- NULL
  if (length(ignored <- intersect(c("startyear", "firstweek", "quarters"), names(cl)))) {
    warning("ignored legacy argument(s): ", paste0(ignored, collapse = ", "))
    cl[ignored] <- NULL
  }
  if (missing(as.one)) cl$as.one <- TRUE  # stsplot_time has different default

  eval.parent(cl)
}
