library("maptools")
library("surveillance")

######################################################################
# Create influenza data for Bavaria and Baden-Wuerttemberg
######################################################################

# read in observed number of cases
flu.counts <- as.matrix(read.table("flu_ByBw.txt"))
namesLK <- substring(colnames(flu.counts),first=2,last=100)
colnames(flu.counts) <- namesLK
# Load population size from table
pop <- as.matrix(read.table("population_2001-12-31_ByBw.txt",header=TRUE)[,c("id","popFrac")])
# Make a matrix containing the population size from the above
# with one row per row in flu.counts 
popM <- matrix(pop[,2],dimnames=list(NULL,pop[,1]),nrow=nrow(flu.counts),ncol=nrow(pop),byrow=TRUE)

# read in adjacency matrix with elements 1 if two regions share a common border
nhood <- as.matrix(read.table("neighourhood_ByBw.txt"))
#map <- readShapePoly("../inst/shapes/districts_BYBW.shp", IDvar = "id")
#Read the shapefile
file <- file.path(.path.package("surveillance"),"shapes","districts_BYBW.shp")
map <- readShapePoly(file, IDvar = "id")


#Create the sts object
fluBB <- new("sts", epoch = 1:nrow(flu.counts),
           observed = flu.counts,
           start = c(2001, 1),
           freq = 52,
           neighbourhood = nhood,
           map = map,
           population = popM
           )

#Spatial plot showing the number of cases in each region for the year 2001
plot(fluBB[year(fluBB) == 2001, ], type= observed ~ 1 | unit , labels = FALSE)
