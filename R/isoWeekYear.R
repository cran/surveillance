######################################################################
# Extract numerical ISO week and year from a Date object
#
# Details:
# This now simply wraps strftime(x, "%V") and strftime(x, "%G"),
# supported on Windows since R 3.1.0. Thus, a handmade implementation
# of isoWeekYear as in surveillance <= 1.16.2 is no longer necessary.
#
# Parameters:
#  Y -- year or a Date/POSIXt object
#  M -- month (only used if Y is the year)
#  D -- day (only used if Y is the year)
#
# Returns:
#  numeric ISO year and week of the date
######################################################################


isoWeekYear <- function(Y, M, D)
{
  if (!inherits(Y, c("Date", "POSIXt")))
    Y <- strptime(paste(Y,M,D,sep="-"),"%Y-%m-%d")

  Wn <- as.numeric(strftime(Y, "%V"))
  Yn <- as.numeric(strftime(Y, "%G"))

  return(list(ISOYear = Yn, ISOWeek = Wn))
}


######################################################################
# An extension of format.Date with additional formatting strings
# - "%Q" / "%OQ" for the quarter (1-4 / I-IV) the month belongs to
# - "%q" days within quarter
# If these formats are not used, base format() is called.
#
# Params:
# x - An object of type Date to be converted.
# format - A character string.
######################################################################

#Small helper function - vectorized gsub, but disregarding names of x
gsub2 <- function(pattern, replacement, x)
{
    len <- length(x)
    mapply(FUN = gsub,
           pattern = rep_len(as.character(pattern), len),
           replacement = rep_len(as.character(replacement), len),
           x = x,
           MoreArgs = list(fixed = TRUE),
           SIMPLIFY = TRUE, USE.NAMES = FALSE)
}

formatDate <- function(x, format)
{
  ##Anything to do?
  if (!grepl("%Q|%OQ|%q", format)) { #nope
    return(format(x,format))
  }

  ##Replicate string
  formatStr <- rep_len(format,length(x))

  ##If days within quarter requested (this is kind of slow)
  if (grepl("%q",format)) {
    ##Loop over vectors of dates
    dateOfQuarter <- sapply(x, function(date) {
      ##Month number in quarter
      modQ <- (as.numeric(format(date,"%m"))-1) %% 3
      dateInMonth <- seq(date,length.out=2,by=paste0("-",modQ," month"))[2]
      ##Move to first of month
      return(dateInMonth - as.numeric(format(dateInMonth,"%d")) + 1)
    })
    dayInQuarter <- as.numeric(x - dateOfQuarter) + 1
    formatStr <- gsub2("%q",as.character(dayInQuarter),formatStr)
  }

  if (grepl("%Q|%OQ",format)) {
    Q <- (as.numeric(format(x,"%m"))-1) %/% 3 + 1 #quarter
    formatStr <- gsub2("%Q",as.character(Q),formatStr)
    formatStr <- gsub2("%OQ",as.roman(Q),formatStr)
  }

  ##The rest of the formatting - works normally as defined by strptime
  res <- character(length(x))
  for (i in 1:length(x))
    res[i] <- format(x[i],formatStr[i])
  return(res)
}
