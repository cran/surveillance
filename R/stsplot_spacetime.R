################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Old implementation of (animated) maps of an sts-object
###
### Copyright (C) 2007-2013 Michael Hoehle
### $Revision: 883 $
### $Date: 2014-04-05 17:26:48 +0200 (Sat, 05 Apr 2014) $
################################################################################


stsplot_spacetime <- function(
    x, type, legend=NULL, opts.col=NULL, labels=TRUE,
    wait.ms=250, cex.lab=0.7, verbose=FALSE, dev.printer=NULL, ...)
{
  #Extract the mappoly
  if (length(x@map) == 0)
    stop("The sts object has an empty map.")
  map <- x@map
  maplim <- list(x=bbox(map)[1,],y=bbox(map)[2,])

  #Check colnames, otherwise no need to continue
  if (is.null(colnames(x@observed)))
    stop("The sts observed slot does not have any colnames to match with the shapefile.")

  #Check for color options
  if (is.null(opts.col)) {
    opts.col <- list(ncolors=100,use.color=TRUE)
  }
  #Check for legend options
  if (is.null(legend)) {
    legend <- list(dx=0.4,dy=0.04,x=maplim$x[1],y=maplim$y[1],once=TRUE)
  }

  #Process dev.printer options
  if (!is.null(dev.printer)) {
    #Device
    if (is.null(dev.printer$device)) dev.printer$device <- png
    #File extension
    if (is.null(dev.printer$extension)) dev.printer$extension <- ".png"
    #Width and height
    if (is.null(dev.printer$width)) dev.printer$width <- 640
    if (is.null(dev.printer$height)) dev.printer$height <- 480
  }
      

  #Extract the data
  o <- x@observed
  alarm <- x@alarm
  
  #Formula is of type "observed ~ 1|unit" (i.e. no time)
  aggregate <- type[[3]][[3]] == "unit"
  if (aggregate) {
    o <- t(as.matrix(apply(o,MARGIN=2,sum)))
    alarm <- t(as.matrix(apply(alarm,MARGIN=2,sum)))>0
  }
  
  #Number of time points
  maxt <- dim(o)[1]

  
  #Get color vector
  gyr <- hcl.colors(ncolors=length(o), use.color=TRUE)
  theCut <- cut(o, length(gyr))
  
  #Cut into specified number of colors
  o.cut <- matrix(as.numeric(theCut),nrow=nrow(o),ncol=ncol(o))
  o.col <- matrix(gyr[o.cut],ncol=ncol(o.cut))
  o.col[is.na(o.col)] <- gray(1)
  dimnames(o.col) <- dimnames(o)

  #Sort the o according to the names in the map
  region.id <- row.names(map)
  o.col.id <- dimnames(o.col)[[2]]

  #Make the columns of o as in the map object
  o.col <- o.col[,pmatch(region.id,o.col.id),drop=FALSE]
  alarm.col <- alarm[,pmatch(region.id,o.col.id),drop=FALSE]

  #Screen processing
  screen.matrix <- matrix(c(0,1,0,1,0,1,0.8,1),2,4,byrow=TRUE)
  split.screen(screen.matrix)

  #Loop over all time slices
  for (t in 1:maxt) {
    #Status information
    if (verbose) {
      cat(paste("Processing slice",t,"of",maxt,"\n"))
    }
    
    #Clean screen (title area)
    screen(n=2)
    par(bg=gray(1))
    erase.screen()
    par(bg="transparent")

    #Plot the map on screen 1
    screen(n=1)
    plot(map,col=o.col[t,],xlab="",ylab="",...)
    #Indicate alarms as shaded overlays
    if (!all(is.na(alarm.col))) {
      #Plotting using density "NA" does not appear to work
      #anymore in the new sp versions
      alarm.col[is.na(alarm.col)] <- 0
      plot(map,dens=alarm.col*15,add=TRUE)
    }
    

    if (labels)
      #getSpPPolygonsLabptSlots is deprecated. Use coordinates method insteas
      text(coordinates(map), labels=as.character(region.id), cex.lab=cex.lab)
  
    if (!aggregate) { title(paste(t,"/",maxt,sep="")) }

    #In case a legend is requested
    if (is.list(legend) && !(legend$once & t>1)  | (t==1)) {
      add.legend(legend, maplim,
                 list(col=gyr, min=min(o), max=max(o), trans=identity))
    }

    #Is writing to files requested?
    if (!is.null(dev.printer)) {
      #Create filename
      fileName <- paste(dev.printer$name,"-",insert.zeroes(t,length=ceiling(log10(maxt))),
                        dev.printer$extension,sep="")
      cat("Creating ",fileName,"\n")
      #Save the current device using dev.print
     dev.print(dev.printer$device, file=fileName, width=dev.printer$width, height=dev.printer$height)
    }
    
    wait(wait.ms) 
  }
  close.screen(all.screens = TRUE)
}



#######################
### auxiliary functions
#######################


### wait a specific amount of milliseconds (via "while" and "proc.time")

wait <- function (wait.ms) # number of milliseconds to wait
{
  #Initialize
  start.time <- proc.time()[3]*1000
  ellapsed <- proc.time()[3]*1000 - start.time

  #Loop as long as required.
  while (ellapsed < wait.ms) {
    ellapsed <- proc.time()[3]*1000 - start.time
  }
}


### add the color key

add.legend <- function(legend, maplim, theColors)
{
  #Preproc
  dy <- diff(maplim$y) * legend$dy
  dx <- diff(maplim$x) * legend$dx
    
  #Add legend -- i.e. a slider
  xlu <- xlo <- legend$x
  xru <- xro <- xlu + dx 
  yru <- ylu <- legend$y
  yro <- ylo <- yru + dy 

  
  step <- (xru - xlu)/length(theColors$col)
  for (i in 0:(length(theColors$col) - 1)) {
    polygon(c(xlo + step * i, xlo + step * (i + 1), 
              xlu + step * (i + 1), xlu + step * i), c(ylo, 
                                                       yro, yru, ylu), col = theColors$col[i + 1], 
            border = theColors$col[i + 1])
  }
  
  
  #Write info about min and max on the slider.
  black <- grey(0)
  lines(c(xlo, xro, xru, xlu, xlo), c(ylo, yro, yru, ylu, ylo), col =   black)

  #Transformation function for data values, e.g., exp or identity
  trans <- theColors$trans

  text(xlu, ylu - 0.5*dy, formatC(trans(theColors$min)), cex = 1, col = black,adj=c(0,1))
  text(xru, yru - 0.5*dy, formatC(trans(theColors$max)), cex = 1, col = black,adj=c(1,1))
}


### Insert leading zeros so integers obtain a fixed length.
### Useful for filenames so they are all of same length (for sorting).

insert.zeroes <- function (x, length=3) # x is an integer
{
  for (i in 1:(length-1)) {
    if (x<10^i) return(paste(paste(rep(0,length-i),collapse=""),x,sep=""))
  }
  #If x has more digits than length then just return x
  return(paste(x))
}
