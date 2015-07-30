buffer.frags <- function(XY, radius, vegetation.base.raster)
#
#Description --
#
#Batch calculates habitat fragmentation indices (‘fragstats’) within a circular buffer zone around sites.
#
#Usage --
#
#source("buffer.frags.R")
#eg.rast <- raster()
#extent(eg.rast) <- c(0,10, 0, 10)
#res(eg.rast) <- c(0.5,0.5)
#eg.rast[] <- NA
#eg.rast[sample(c(1:ncell(eg.rast)), 100)] <- 1
#plot(eg.rast)
#eg.coords <- data.frame(Longitude=c(3, 7), Latitude=c(4, 8))
#row.names(eg.coords) <- c("SiteA", "SiteB")
#eg.coords
#points(eg.coords, pch=20, cex=2)
#frags <- buffer.frags(eg.coords, radius=150000, vegetation.base.raster=eg.rast)
#
#frags <- do.call(rbind, frags) #TO CONDENSE INTO SINGLE DATA.FRAME
#frags #NOTE:class =100 in 'frags' is the default for NAs in the input raster (i.e. no habitat/vegetation).
#
#Arguments --
#
#XY: vector referring to  the x/longitude and y/latitude coordinates of a set of sites. If the coordinates are in columns of a data.frame, then the subsetted data.frame should be given: e.g. XY=df[,"Longitude", "Latitude"]. If the coordinates are within a SpatialPointsDataFrame, the SPDF can be given for XY (function will locate the coordinates)
#
#radius: the radius of the buffer zone around each site within which fragstats are calculated, in units of metres if raster is in longlat.
#
#vegetation.base.raster: a rasterLayer object containing integer classes for the presence of vegetation/habitat of 1 or more kinds.
#
#Details --
#
#This function is a wrapper for the ClassStat function from package SDMTools, which calculates fragstats for a given raster representing a vegetation/habitat matrix. The ClassStat function is applied to a user defined circular area surrounding a site/coordinate, and these are batch processed for multiple sites. The circular buffer is generated through the function raster::buffer, which can be slow for high resolution rasters or large areas. 
#
#Value --
#
#Returns a list of data frames, one for each focal location, with fragstat indices as columns and vegetation/habitat classes (read from the input raster) as rows. See example for collapsing the list to a single data.frame object. By default, NA values within the buffer (generally representing non-habitat or non-vegetation) are returned with the value 100 for 'class'. 'Class' with otherwise be the same as the categories in the habitat raster (See ?ClassStat)
#
#Required packages --
#
#raster, SDMTools
#
#Authors --
#
#Greg R. Guerin
#
#References --
#
#Wrapper for SDMTools version of fragstats, which is based on statistics calculated by fragstats, see http://www.umass.edu/landeco/research/fragstats/fragstats.html.
#
#Jeremy VanDerWal, Lorena Falconi, Stephanie Januchowski, Luke Shoo and Collin Storlie (2014). SDMTools: Species Distribution Modelling Tools: Tools for processing data associated with species distribution modelling exercises. R package version 1.1-221. http://CRAN.R-project.org/package=SDMTools
#
#
#License --
#
#GNU GPL-3
#
#Version --
#
#1.0
#
{
	if(class(XY) == "SpatialPointsDataFrame") {
		LONG <- coordinates(XY)[,1]
		LAT <- coordinates(XY)[,2]
	}
	
	
	if(class(XY) == "data.frame") {
		LONG <- XY[,1]
		LAT <- XY[,2]
	}
	
	if(class(XY) != "SpatialPointsDataFrame" & class(XY) != "data.frame") {
		stop("XY must be a data.frame with x/longitude and y/latitude columns or a SpatialPointsDataFrame")
	}
	
	
	if(class(radius) != "numeric") {
		stop("'radius' should be a number indicating the radius of the buffer.")
	}
	
	if(class(vegetation.base.raster) != "RasterLayer") {
		stop("'vegetation.base.raster' should be a RasterLayer representing vegetation or habitat presence.")
	}
	
	require(raster)
	require(SDMTools)
	
	n <- 0
	dat <- list()
	for(i in 1:length(LONG)) {
		n <- n + 1
		temp.rast <- vegetation.base.raster #copy
		focal.cell <- cellFromXY(temp.rast, c(LONG[n], LAT[n])) #find raster cell associated with the site/coordinates
		temp.rast[] <- NA
		temp.rast[focal.cell] <- 1
		temp.rast <- raster::buffer(temp.rast, width=radius)
		outer.cells <- which(is.na(temp.rast[])) #vector of cells outside the buffer zone
		temp.rast <- vegetation.base.raster #reset
		temp.rast[is.na(temp.rast)] <- 100
		temp.rast[outer.cells] <- NA
		dat[[n]] <- ClassStat(temp.rast)
		}
		try(names(dat) <- row.names(XY))
		return(dat)
	}
