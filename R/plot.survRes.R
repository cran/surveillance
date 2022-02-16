## legacy code removed in surveillance 1.20.0
## now plotting via stsplot_time()
plot.survRes <- function(x, method = x$control$name, disease = x$control$data,
                         xaxis.years = TRUE, startyear = 2001, firstweek = 1,
                         same.scale = TRUE, ...,
                         main = paste0("Analysis of ", disease, " using ", method))
{
  stopifnot(is.vector(x$control$range, mode = "numeric"))
  stsObjRange <- disProg2sts(x[["disProgObj"]])[x$control$range,]
  stsObjRange@alarm[] <- x[["alarm"]]
  stsObjRange@upperbound[] <- x[["upperbound"]]

  if (length(ignored <- c("startyear", "firstweek")[c(!missing(startyear), !missing(firstweek))]))
    warning("ignored legacy argument(s): ", paste0(ignored, collapse = ", "))

  ## see plot(sts.cdc) in vignette("surveillance") or demo("cost")
  hookFunc <- if (is.null(x$aggr)) function() NULL else {
    stsObjRange@control$aggr <- x$aggr
    upperboundx <- NULL  # make codetools::checkUsage() happy
    function() points(upperboundx, x@control$aggr)
  }

  if (xaxis.years) {
      stsplot_time(stsObjRange, same.scale = same.scale, ..., main = main,
                   hookFunc = hookFunc)
  } else {
      stsplot_time(stsObjRange, same.scale = same.scale, ..., main = main,
                   hookFunc = hookFunc, xaxis.labelFormat = NULL)
  }
}
