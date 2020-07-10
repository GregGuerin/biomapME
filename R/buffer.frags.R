buffer.frags <- function(XY, radius, vegetation.base.raster, plot=TRUE) {
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


  if(1 %in% grep("longlat", projection(vegetation.base.raster))) {longlat=TRUE} else { 
  longlat=FALSE}
  buffers <- dismo:::.generateCircles(XY, d=radius, lonlat=longlat)
  vegetation.base.raster[is.na(vegetation.base.raster)] <- 100
  
  
  if(plot){
    plot(vegetation.base.raster)
    points(XY, pch=20)
    plot(buffers, add=TRUE)
  } #cls if plot
  
	n <- 0
	dat <- list()
	for(i in 1:length(LONG)) {
		n <- n + 1
		# temp.rast <- vegetation.base.raster #copy
		# focal.cell <- cellFromXY(temp.rast, c(LONG[n], LAT[n])) #find raster cell associated with the site/coordinates
		# temp.rast[] <- NA
		# temp.rast[focal.cell] <- 1
		# temp.rast <- raster::buffer(temp.rast, width=radius)
		# outer.cells <- which(is.na(temp.rast[])) #vector of cells outside the buffer zone
		# temp.rast <- vegetation.base.raster #reset
		# temp.rast[is.na(temp.rast)] <- 100
		# temp.rast[outer.cells] <- NA
		# plot(temp.rast)
		temp.rast <- mask(vegetation.base.raster, buffers[n])
		#if(plot) {dev.new()
		#plot(temp.rast)}
		dat[[n]] <- ClassStat(temp.rast)
		}
		try(names(dat) <- row.names(XY))
		return(dat)
	}
