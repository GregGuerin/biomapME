range.metrics <- function(species_records, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", coord.type="longlat", weight.type="cell", geo.calc="max.dist", outlier_pct=95, verbose=TRUE, frame.raster, deg.resolution=c(0.25,0.25), extent.vector, plot.out=TRUE) {

# 	require(adehabitat)

	if(outlier_pct > 99 | outlier_pct < 1) {
		stop("Outlier_pct should be a percentage")
	} #cls outlier stop...

	if(class(species_records) != "data.frame") {
		stop("Species data must be in a data.frame")
	} #cls species_records class check




		species_records <- species_records[,c(species, longitude, latitude)]
		colnames(species_records) <- c("SPECIES", "LONGITUDE", "LATITUDE")


	if(!("SPECIES" %in% colnames(species_records))) {stop("Cannot locate species data")}
	if(!("LATITUDE" %in% colnames(species_records))) {stop("Cannot locate latitude/y data")}
	if(!("LONGITUDE" %in% colnames(species_records))) {stop("Cannot locate longitude/x data")}
	if(any(is.na(species_records$LONGITUDE))) {
		species_records <- species_records[-which(is.na(species_records$LONGITUDE)),]
	} #cls NA longitude
	if(any(is.na(species_records$LATITUDE))) {
		species_records <- species_records[-which(is.na(species_records$LATITUDE)),]
	} #cls NA latitude

	coordinates(species_records) <- c("LONGITUDE", "LATITUDE")

	#####
	if(plot.out == TRUE) {
		if(coord.type == "longlat") {
		  lon.ext <- extent(species_records)[1:2]
		  lat.ext <- extent(species_records)[3:4]
		  map("world", fill=TRUE, col="gray50", bg="lightblue", ylim=lat.ext, xlim=lon.ext, mar=c(0,0,0,0))
		  points(species_records, col=species_records$SPECIES, pch=16)
		} #cls if longlat
	  if(coord.type == "custom") {
	    plot(species_records, col=species_records$SPECIES, pch=16)
	  } #cls if custom
		} #cls if(plot.out)
	######


	if(weight.type=="cell") {



	if(missing(frame.raster)) {
		frame.raster <- raster()
		if(missing(extent.vector)) {
			extent(frame.raster)@xmin <- floor(extent(species_records)@xmin)
			extent(frame.raster)@ymin <- floor(extent(species_records)@ymin)
			extent(frame.raster)@xmax <- ceiling(extent(species_records)@xmax)
			extent(frame.raster)@ymax <- ceiling(extent(species_records)@ymax)
			}
		if(!(missing(extent.vector))) {
			extent(frame.raster) <- extent.vector
		}
		res(frame.raster) <- deg.resolution
		cat("Generating frame raster at ", deg.resolution, " resolution and extent defined by: ", extent(frame.raster)@xmin, extent(frame.raster)@xmax, extent(frame.raster)@ymin, extent(frame.raster)@ymax,"\n")
	} #cls if(missing(frame.raster))


	if(!(extent(species_records)@xmin >= extent(frame.raster)@xmin & extent(species_records)@xmax <= extent(frame.raster)@xmax & extent(species_records)@ymin >= extent(frame.raster)@ymin & extent(species_records)@ymax <= extent(frame.raster)@ymax)) {
		cat("Some point locations lie outside the frame raster -- trimming these records", "\n")
		species_record_COORDS <- as.data.frame(coordinates(species_records))
		species_records <- species_records[-which(species_record_COORDS$LONGITUDE < extent(frame.raster)@xmin | species_record_COORDS$LONGITUDE > extent(frame.raster)@xmax | species_record_COORDS$LATITUDE < extent(frame.raster)@ymin | species_record_COORDS$LATITUDE > extent(frame.raster)@ymax),]
		}

	frame.raster[] <- NA



	cell.ranges <- function(x) { #where x is species_records (a SPDF object)
		n <- 0
		v <- rep(0, length(unique(x$SPECIES)))
		for (i in unique(x$SPECIES)) {
			n <- n + 1
			temp <- x[which(x$SPECIES == i),]
			occupied.cells <- cellFromXY(frame.raster, temp)
			frame.raster[occupied.cells] <- 1
			numberOFcells <- sum(frame.raster[!is.na(frame.raster[])])
			v[n] <- numberOFcells
			if(plot.out == TRUE) {
				dev.new()
				plot(frame.raster, main=i, breaks=c(0.5,1.5), col=topo.colors(length(unique(species_records$SPECIES)))[n])
			} #cls if(plot.out)
			names(frame.raster) <-  i
			if(n == 1) {
				maps <- stack(frame.raster)
			}#cls if(n ==1)
			if(n > 1) {
				maps <- stack(maps, frame.raster)
			} #cls if(n >1)
			frame.raster[] <- NA #reset
		} #cls each species loop
		names(v) <- unique(x$SPECIES)
		names(maps) <- unique(x$SPECIES)
		outputs <- list(ranges=v, rasters=maps)
		return(outputs)
	} #cls cell.ranges function


	ranges <- cell.ranges(species_records)



	} #cls if(weight.type)=="cell"








		if(weight.type=="geo") {


			CalcDists <- function(longlats) { #modified from CalcDists.R, see https://gist.githubusercontent.com/sckott/931445/raw/9db1d432b2308a8861f6425f38aaabbce44eb994/CalcDists.R
				name <- list(rownames(longlats), rownames(longlats))
				n <- nrow(longlats)
				z <- matrix(0, n, n, dimnames = name)
				for (i in 1:n) {
					for (j in 1:n) {
						if(coord.type=="longlat") {
							z[i, j] <- distCosine(c(longlats[j, 1], longlats[j, 2]), c(longlats[i, 1], longlats[i, 2]))
						} #cls if(missing(XY))
						if(coord.type=="custom") {
							z[i, j] <- sqrt(sum((c(longlats[j, 1], longlats[j, 2]) - c(longlats[i, 1], longlats[i, 2])) ^ 2))
						} #cls if(!(missing(XY)))
						} #cls for (j)
					} #cls for (i)
				z <- as.dist(z)
				return(z)
				} #cls CalcDists





				spp_ranges <- function(x) {
					n <- 0
					v <- rep(0, length(unique(x$SPECIES)))
					x$SPECIES <- sub(pattern = " ", replacement = ".", x = x$SPECIES, fixed=TRUE)
					x$SPECIES <- sub(pattern = "-", replacement = ".", x = x$SPECIES, fixed=TRUE)
					for (i in unique(x$SPECIES)) {
						n <- n + 1
						temp <- x[which(x$SPECIES == i),]
						colnames(temp) <- colnames(x) #?necessary

								if(geo.calc == "max.dist") {
									if(nrow(temp) < 2) {
										v[n] <- 1
										warning("Only one record, returning a default range size of 1 for ", i)
									} else {
											if(nrow(temp) < 5) {
											  temp2 <- data.frame(LONGITUDE=temp$LONGITUDE, LATITUDE=temp$LATITUDE)
											  v[n] <- max(CalcDists(temp2))
											} #cls if(nrow(temp) < 5)...
											if(nrow(temp) > 4) {
												spp_i_range_polygon <- try(mcp(temp, percent=outlier_pct))
												if(class(spp_i_range_polygon)[1] == "try-error") {
												  temp2 <- data.frame(LONGITUDE=temp$LONGITUDE, LATITUDE=temp$LATITUDE)
												  v[n] <- max(CalcDists(temp2))
												} #cls if(class...
												if(!class(spp_i_range_polygon)[1] == "try-error") {
													v[n] <- max(CalcDists(as.data.frame(spp_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords)))
												} #cls if(!(class...
											} #cls if(nrow(temp) > 4)...
										} #cls else...
								} #cls if(geo.calc == "max.dist")...

								if(geo.calc == "polygon") {
									if(nrow(temp) < 5) {
										if(coord.type=="longlat") {
											polygon_area <- try(areaPolygon(temp))
											} #cls if(coord.type="longlat")
										if(coord.type=="custom") {
											polygon_area <- try(abs(polyarea(temp[,"LONGITUDE"], temp[,"LATITUDE"])))
										} #cls if(coord.type="custom")
										if(class(polygon_area) == "try-error") {
											v[n] <- 1
											warning("Cannot compute polygon, returning 1 as the default range area for ", i)
										} #cls if(class(polygon_area)...
										if(class(polygon_area) == "numeric") {
											if(polygon_area == 0) {
												v[n] <- 1
												warning("Cannot compute polygon, returning 1 as the default range area for ", i)
											} #cls is(polygon_area == 0)
											if(polygon_area != 0) {
												v[n] <- polygon_area
											} #cls if(polygon_area != 0)

										} #cls if(class(polyon_area) == "numeric...
									} #cls if(nrow(temp) < 5)...
									if(nrow(temp) > 4) {
										spp_i_range_polygon <- try(mcp(temp, percent=outlier_pct))
										if(class(spp_i_range_polygon)[1] == "try-error") {
											v[n] <- 1
											warning("Unable to compute a polygon, returning 1 as the range area for ", i)
										} #cls if(class(spp...
										if(class(spp_i_range_polygon)[1] != "try-error") {
											if(coord.type=="longlat") {
												polygon_area <- try(areaPolygon(as.data.frame(spp_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords)))
											} #cls if(coord.type="longlat")
											if(coord.type=="custom") {
												polygon_area <- try(abs(polyarea(spp_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords[,"LONGITUDE"], spp_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords[,"LATITUDE"])))
											} #cls if(coord.type="custom")
											if(class(polygon_area) == "try-error") {
												v[n] <- 1
												warning("Cannot compute polygon area, returning 1 as the range area for ", i)
											} #cls if(class(polygon_area)...
											if(class(polygon_area) == "numeric") {
												v[n] <- polygon_area
												if(plot.out == TRUE) {
												  plot(spp_i_range_polygon, main=paste(polygon_area))
												} #cls if plot.out...
											} #cls if(class(poygon_area == numeric...
										} #cls if(!class(spp...
									} #cls if(nrow(temp) > 4)...
								} #cls if geo.calc = polygon


						if(geo.calc == "LONG") {
						  v[n] <- max(temp$LONGITUDE)-min(temp$LONGITUDE)
						} #cls if(geo.calc =="LONG")...

						if(geo.calc == "LAT") {
						  v[n] <- max(temp$LATITUDE)-min(temp$LATITUDE)
						} #cls if(geo.calc =="LONG")...


						if(verbose) cat(n, "species complete:", i, v[n], "\n")

						} #cls for (i in colnames...

					names(v) <- unique(x$SPECIES)
					return(v)
					} #cls spp_ranges function

				ranges <- spp_ranges(species_records)







		} #close if(weight.type = geo...



	##collate outputs

		if(weight.type == "geo") {
		  outputs <- list(weights = ranges)
		  return(outputs)
		} #cls if "geo"...

	  if(weight.type == "cell") {
	    outputs <- list(weights = ranges$ranges, rasters = ranges$rasters)
	    return(outputs)
	  } #cls if "cell"...

} #cls function


