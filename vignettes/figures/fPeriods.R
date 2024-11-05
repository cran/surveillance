#!/usr/bin/env Rscript
# Create plots for fig:fPlot in vignette("monitoringCounts")
# Author: Ma"elle Salmon
# Code removed from vignette to drop ggplot2 dependency
library(ggplot2)

# for rectangles
widthRectangles <- 10
# dimensions for the ticks
heightTick <- 4
xTicks <- c(15,67,119)
yTicksStart <- rep(0,3)
yTicksEnd <- rep(0,3)
yTicksEnd2 <- rep(-5,3)
textTicks <- c("t-2*p","t-p","t[0]")
xBigTicks <- c(xTicks[1:2]-widthRectangles/2,xTicks[1:2]+widthRectangles/2,xTicks[3]-widthRectangles/2,xTicks[3])
yTicksBigEnd <- rep(0,6)
yTicksBigStart <- rep(heightTick,6)
# to draw the horizontal line
vectorDates <- rep(0,150)
dates <- seq(1:150)
data <- data.frame(dates,vectorDates)
xPeriods <- c(15,67,117,15+26,67+26)

p <- ggplot() +
# white
  theme_void() +
  geom_segment(aes(x = 0, y = -20, xend = 200, yend = 10), linewidth = 2,
               arrow = arrow(length = unit(0.5, "cm")), colour ='white')   +
# time arrow
  geom_segment(aes(x = 0, y = 0, xend = 150, yend = 0), linewidth = 1,
               arrow = arrow(length = unit(0.5, "cm")))   +
# ticks
  geom_segment(aes(x = xTicks, y = yTicksEnd2, xend = xTicks, yend = yTicksStart ),
               arrow = arrow(length = unit(0.3, "cm")), linewidth = 1) +
# big ticks
  geom_segment(aes(x = xBigTicks, y = yTicksBigStart, xend = xBigTicks, yend = yTicksBigEnd*2),
               linewidth = 1) +
# time label
  annotate("text", label = "Time", x = 170, y = 0, size = 8, colour = "black",
           family="serif") +
# ticks labels
  annotate('text',label=c("t[0]-2 %.% freq","t[0]-freq","t[0]"),x = xTicks,
           y = yTicksEnd - 10, size = 8,family="serif",parse=T)

## noPeriods = 2
pdf("fPeriods2.pdf", width = 7, height = 3, colormodel = "gray")
p +
# periods labels
annotate('text',label=c("A","A","A","B","B"),x = xPeriods,
         y = rep(6,5), size = 8,family="serif",parse=T)
dev.off()

## noPeriods = 3
yTicksBigEnd2 <- rep(0,4)
yTicksBigStart2 <- rep(heightTick,4)
newX <- c(xTicks[1:2]+widthRectangles/2+52-widthRectangles,xTicks[1:2]+52/2)
xPeriods <- c(15,67,117,15+16,67+16,15+35,67+35)
pdf("fPeriods3.pdf", width = 7, height = 3, colormodel = "gray")
p +
  geom_segment(aes(x = newX, y = yTicksBigStart2, xend = newX, yend = yTicksBigEnd2),
               linewidth = 1) +
# periods labels
annotate('text',label=c("A","A","A","B","B","C","C"),x = xPeriods,
         y = rep(6,7), size = 8,family="serif",parse=T)
dev.off()
