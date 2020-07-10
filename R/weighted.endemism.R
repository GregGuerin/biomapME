weighted.endemism <- function(species_records, records="single", site.coords, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", frame.raster, deg.resolution=c(0.25,0.25), extent.vector, type="weighted", plot.raster=TRUE, own.weights, weight.type="cell", geo.type="cell", geo.calc="max.dist", outlier_pct=100, verbose=TRUE, own.grid.matrix) {

	if(outlier_pct > 100 | outlier_pct < 1) {
		stop("Outlier_pct should be a percentage")
	}

	if(class(species_records) != "data.frame") {
		stop("Species data must be in a data.frame")
	}

	if(records == "site") {
		species_records <- convert.site.data(species_records, site.coords)
	}

	if(records == "single") {
		species_records <- species_records[,c(species, longitude, latitude)]
		colnames(species_records) <- c("SPECIES", "LONGITUDE", "LATITUDE")
	}


	if(!("SPECIES" %in% colnames(species_records))) {stop("Cannot locate species data")}
	if(!("LATITUDE" %in% colnames(species_records))) {stop("Cannot locate latitude data")}
	if(!("LONGITUDE" %in% colnames(species_records))) {stop("Cannot locate longitude data")}
	if(any(is.na(species_records$LONGITUDE))) {
		species_records <- species_records[-which(is.na(species_records$LONGITUDE)),]
	}
	if(any(is.na(species_records$LATITUDE))) {
		species_records <- species_records[-which(is.na(species_records$LATITUDE)),]
	}

	coordinates(species_records) <- c("LONGITUDE", "LATITUDE")

	if(missing(frame.raster)) {
		if(!(missing(own.grid.matrix))) {stop("You must supply a frame.raster with cells that match own.grid.matrix")}
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
	}


	if(!(extent(species_records)@xmin >= extent(frame.raster)@xmin & extent(species_records)@xmax <= extent(frame.raster)@xmax & extent(species_records)@ymin >= extent(frame.raster)@ymin & extent(species_records)@ymax <= extent(frame.raster)@ymax)) {
		cat("Some point locations lie outside the frame raster -- trimming these records", "\n")
		species_record_COORDS <- as.data.frame(coordinates(species_records))
		species_records <- species_records[-which(species_record_COORDS$LONGITUDE < extent(frame.raster)@xmin | species_record_COORDS$LONGITUDE > extent(frame.raster)@xmax | species_record_COORDS$LATITUDE < extent(frame.raster)@ymin | species_record_COORDS$LATITUDE > extent(frame.raster)@ymax),]
		}


	if(!(missing(own.grid.matrix))) {
		cat("Reading user-defined gridded occurrence matrix", "\n")
		if(class(own.grid.matrix) != "data.frame") {stop("Supplied own.grid.matrix must be a data.frame")}
		cell_occur_matrix <- own.grid.matrix
	} #clse if(!(missing(own.grid...

	if(missing(own.grid.matrix)) {
	  cell_occur_matrix <- map.pa.matrix(as.data.frame(species_records), frame.raster=frame.raster)$grid.matrix
	} #cls if(missing(own.grid...




	if(!(missing("own.weights"))) {
		if(class(own.weights) != "numeric") {stop("Supplied species weights must be numeric")}
		if(!(all(names(own.weights) %in% colnames(cell_occur_matrix)) & all(colnames(cell_occur_matrix) %in% names(own.weights)))) {
		  stop("Species names for supplied weights are not a complete match for species in supplied records")
		  } #true if all in it both ways
		if(length(own.weights) != length(colnames(cell_occur_matrix))) {
		  stop("Supplied weights vector is different length than number of species - must supply a weighting for each species")
		  }

		cat("Calculating user supplied weights", "\n")
		inv_rang_cell_occur_mat <- cell_occur_matrix
		for (i in 1:ncol(inv_rang_cell_occur_mat)) {inv_rang_cell_occur_mat[,i] <- inv_rang_cell_occur_mat[,i]/own.weights[which(names(own.weights) == colnames(inv_rang_cell_occur_mat)[i])]}
		ranges <- own.weights
	} #cls if(!(missing...



	if(missing("own.weights")) {


		if(weight.type=="cell") {
			cat("Calculating cell-based range weights", "\n")
			inv_rang_cell_occur_mat <- apply(cell_occur_matrix, 2, function(x) {x/sum(x)})
			ranges <- colSums(cell_occur_matrix)
		} #cls if(weight.type = cell...

		if(weight.type=="richness") {
			inv_rang_cell_occur_mat <- cell_occur_matrix
			ranges <- apply(cell_occur_matrix, 2, function(x) {x = 1})
		} #cls if(w.t = richness...

		if(weight.type=="geo") {

		  #####################################
			CalcDists <- function(longlats) { #modified from CalcDists.R, see https://gist.githubusercontent.com/sckott/931445/raw/9db1d432b2308a8861f6425f38aaabbce44eb994/CalcDists.R
				name <- list(rownames(longlats), rownames(longlats))
				n <- nrow(longlats)
				z <- matrix(0, n, n, dimnames = name)
				for (i in 1:n) {
					for (j in 1:n) z[i, j] <- distCosine(c(longlats[j, 1], longlats[j, 2]), c(longlats[i, 1], longlats[i, 2]))
					}
				z <- as.dist(z)
				return(z)
				} #cls CalcDists
      #####################################

			if(geo.type == "cell") {

				cell_centroids <- as.data.frame(coordinates(frame.raster))
				colnames(cell_centroids) <- c("LONGITUDE", "LATITUDE")
				cell_centroids$cells <- row.names(cell_centroids)

				species_ranges <- function(x) {
					v <- rep(0, ncol(x))
					for (i in 1:ncol(x)) {
						temp <- as.data.frame(x[,i])
						colnames(temp) <- "species_i"
						row.names(temp) <- row.names(x)
						temp$delete <- temp$species_i
						if(any(temp$species_i == 0)) {
							temp <- temp[-which(temp$species_i == 0),]
							} #cls if(any(temp...
							temp$cells <- row.names(temp)
							temp <- merge(temp, cell_centroids, by="cells")[,c("LONGITUDE", "LATITUDE")]
							coordinates(temp) <- c("LONGITUDE", "LATITUDE") #xy of occupied cell centroids

							if(geo.calc == "max.dist") {
								cell_dimensions <- (mean(values(raster::area(frame.raster)))*1000000)^(0.5) #area in km^2 so converting to m^2 to match areaPolygon, then find square root to get back to 1-dimension in m (it is ~110 km per)
								if(nrow(as.data.frame(temp)) < 2) {
									v[i] <- cell_dimensions
									warning("Only one record, returning the span of single grid cell for ", colnames(x[i]))
									} else {
										if(nrow(as.data.frame(temp)) < 5) {
											v[i] <- max(CalcDists(as.data.frame(temp)))
											if(v[i] < cell_dimensions) {
												v[i] <- cell_dimensions
												warning("Range is less than spatial grain of frame.raster, returning span of single grid cell for ", colnames(x[i]))
											} #cls if(v[i] < cell_dimensions)...
										} #cls if nrow(temp) < 5...
										if(nrow(as.data.frame(temp)) > 4) {
											spp_i_range_polygon <- try(mcp(temp, percent=outlier_pct))
											if(class(spp_i_range_polygon)[1] == "try-error") {
												v[i] <- max(CalcDists(as.data.frame(temp)))
											} #cls if(class(spp...
											if(!class(spp_i_range_polygon)[1] == "try-error") {
												v[i] <- max(CalcDists(as.data.frame(spp_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords)))
											}	 #cls if(!class(spp...
											if(v[i] < cell_dimensions) {
												v[i] <- cell_dimensions
												warning("Range is less than spatial grain of frame.raster, returning span of single grid cell for ", colnames(x[i]))
												} #cls if(v[i] <...
										}	#cls if(nrow(temp) > 4)...
									} #cls else...
							} #cls if geo.calc = max.dist...

							if(geo.calc == "polygon") {
								cell_size <- mean(values(raster::area(frame.raster)))*1000000 #area in km^2 so convert to m^2 to match areaPolygon
								if(nrow(as.data.frame(temp)) < 5) {
									polygon_area <- try(areaPolygon(as.data.frame(temp)))
									if(class(polygon_area) == "try-error") {
										v[i] <- cell_size
										warning("Less than minimum number of required points to compute polygon, returning approximate area of single grid cell for ", colnames(x[i]))
									} #cls if(class(polyon...
									if(class(polygon_area) == "numeric") {
										if(polygon_area < cell_size) {
											v[i] <- cell_size
											warning("Polygon area is less than spatial grain of frame.raster, returning approximate area of a single grid cell for ", colnames(x[i]))
										} #cls if(polygon_area <...
										if(polygon_area >= cell_size) {
											v[i] <- polygon_area
										} #CLS if(class(polygon..
									} #cls if(class(polygon...
								} #cls if(nrow(temp) < 5)...
								if(nrow(as.data.frame(temp)) > 4) {
									spp_i_range_polygon <- try(mcp(temp, percent=outlier_pct))
									if(class(spp_i_range_polygon)[1] == "try-error") {
										v[i] <- cell_size
										warning("Unable to compute polygon, returning approximate area of single grid cell for ", colnames(x[i]))
									} #cls if(class(spp...
									if(class(spp_i_range_polygon)[1] != "try-error") {
										polygon_area <- try(areaPolygon(as.data.frame(spp_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords)))
										if(class(polygon_area) != "numeric") {
											v[i] <- cell_size
											warning("Cannot compute polygon area, returning approximate area of single grid cell for ", colnames(x[i]))
										} #cls if(class(polyon...
										if(class(polygon_area) == "numeric") {
											if(polygon_area < cell_size) {
												v[i] <- cell_size
												warning("Polygon area is less than spatial grain of frame.raster, returning approximate area of single grid cell for ", colnames(x[i]))
											} #cls if(polygon_area <...
											if(polygon_area >= cell_size) {
												v[i] <- polygon_area
											} #cls if(polygon_area >=...
										} #cls if(class(polygon...
									} #cls if(!class(spp...
								} #cls if(nrow(temp) > 4)...
							} #cls if(geo.calc = polygon)

					if(verbose) cat("spp complete:", i, v[i], colnames(x[i]), "\n")

				} #cls for (i in...

				names(v) <- colnames(cell_occur_matrix)
				return(v)
				} #cls species_ranges function

				ranges <- species_ranges(cell_occur_matrix)

			} #cls if(geo.type = cell...


			if(geo.type == "point") {

				spp_ranges <- function(x) {
					n <- 0
					v <- rep(0, length(colnames(cell_occur_matrix)))
					x$SPECIES <- sub(pattern = " ", replacement = ".", x = x$SPECIES, fixed=TRUE)
					x$SPECIES <- sub(pattern = "-", replacement = ".", x = x$SPECIES, fixed=TRUE)
					for (i in colnames(cell_occur_matrix)) {
						n <- n + 1
						temp <- x[which(x$SPECIES == i),]


								if(geo.calc == "max.dist") {
									cell_dimensions <- (mean(values(raster::area(frame.raster)))*1000000)^(0.5) #area in km^2 so convert to m^2 to match areaPolygon, then find square root to get back to 1-dimension
									if(nrow(as.data.frame(temp)) < 2) {
										v[n] <- cell_dimensions
										warning("Only one record, returning the size of a single grid cell (1-dimension) for ", i)
									} else {
											if(nrow(as.data.frame(temp)) < 5) {
												v[n] <- max(CalcDists(as.data.frame(temp)[,c("LONGITUDE", "LATITUDE")]))
												if(v[n] < cell_dimensions) {
													v[n] <- cell_dimensions
													warning("Range span less than spatial grain of frame.raster, returning span of single grid cell for ", i)
												} #cls if(v[n] < cell_dimensions...
											} #cls if(nrow(temp) < 5)...
											if(nrow(as.data.frame(temp)) > 4) {
												spp_i_range_polygon <- try(mcp(temp, percent=outlier_pct))
												if(class(spp_i_range_polygon)[1] == "try-error") {
													v[n] <- max(CalcDists(as.data.frame(temp)[,c("LONGITUDE", "LATITUDE")]))
												} #cls if(class...
												if(!class(spp_i_range_polygon)[1] == "try-error") {
													v[n] <- max(CalcDists(as.data.frame(spp_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords)))
												} #cls if(!(class...
												if(v[n] < cell_dimensions) {
													v[n] <- cell_dimensions
													warning("Range span less than spatial grain of frame.raster, returning span of single grid cell for ", i)
												} #cls if(v[n] <...
											} #cls if(nrow(temp) > 4)...
										} #cls else...
								} #cls if(geo.calc == "max.dist")... [three open now...]

								if(geo.calc == "polygon") {
									cell_size <- mean(values(raster::area(frame.raster)))*1000000 #returns size in m^2 (convert from km^2)
									if(nrow(as.data.frame(temp)) < 5) {
										polygon_area <- try(areaPolygon(as.data.frame(temp)[,c("LONGITUDE", "LATITUDE")]))
										if(class(polygon_area) == "try-error") {
											v[n] <- cell_size
											warning("Cannot compute polygon, returning approximate area of single grid cell for ", i)
										} #cls if(class(polygon_area)...
										if(class(polygon_area) == "numeric") {
											if(polygon_area < cell_size) {
												v[n] <- cell_size
												warning("Polygon area less than spatial grain of frame.raster, returning approximate area of single grid cell for ", i)
											} #cls if(polygon_area <...
											if(polygon_area >= cell_size) {
												v[n] <- polygon_area
											} #cls if (polygon_area >=...
										} #cls if(class(polyon_area) == "numeric...
									} #cls if(nrow(temp) < 5)...
									if(nrow(as.data.frame(temp)) > 4) {
										spp_i_range_polygon <- try(mcp(temp, percent=outlier_pct))
										if(class(spp_i_range_polygon)[1] == "try-error") {
											v[n] <- cell_size
											warning("Unable to compute a polygon, returning approximate area of single grid cell for ", i)
										} #cls if(class(spp...
										if(class(spp_i_range_polygon)[1] != "try-error") {
											polygon_area <- try(areaPolygon(as.data.frame(spp_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords)))
											if(class(polygon_area) == "try-error") {
												v[n] <- cell_size
												warning("Cannot compute polygon area, returning approximate area of single grid cell for ", i)
											} #cls if(class(polygon_area)...
											if(class(polygon_area) == "numeric") {
												v[n] <- polygon_area
												if(v[n] < cell_size) {
													v[n] <- cell_size
													warning("Polygon area is less than spatial grain of frame.raster, returning approximate area of single grid cell for ", i)
												} #cls if(v[n] < cell_size)...
											} #cls if(class(poygon_area == numeric...
										} #cls if(!class(spp...
									} #cls if(nrow(temp) > 4)...
								} #cls if geo.calc = polygon

						if(verbose) cat(n, "species complete:", i, v[n], "\n")

						} #cls for (i in colnames...

					names(v) <- colnames(cell_occur_matrix)
					return(v)
					} #cls spp_ranges function

				ranges <- spp_ranges(species_records)

			} #cls if geo.type = point...



			cat("Calculating geographic range weights", "\n")
			inv_rang_cell_occur_mat <- cell_occur_matrix
			for(i in 1:ncol(inv_rang_cell_occur_mat)) {
				inv_rang_cell_occur_mat[,i] <- inv_rang_cell_occur_mat[,i]/ranges[which(names(ranges) == colnames(inv_rang_cell_occur_mat)[i])]
			} #cls for(i in 1:ncol...


		} #close if(weight.type = geo...
	} #close if(missing(own.weights...

	cat("Calculating weighted endemism", "\n")
	rawEndemism <- rowSums(inv_rang_cell_occur_mat)


	if(type=="weighted") {
		WE_raster <- frame.raster
		WE_raster[] <- NA
		WE_raster[as.numeric(names(rawEndemism))] <- rawEndemism
		if(plot.raster==TRUE) {plot(WE_raster, main="Weighted Endemism")}
		outputs <- list(WE = rawEndemism, WE_raster = WE_raster, weights = ranges, grid.matrix = cell_occur_matrix)
		return(outputs)
	} #cls if(type = w...

	if(type=="corrected"){
		richness4endemism <- rowSums(cell_occur_matrix)
		corrected.WE <- rawEndemism/richness4endemism
		corrected.WE_raster <- frame.raster
		corrected.WE_raster[] <- NA
		corrected.WE_raster[as.numeric(names(corrected.WE))] <- corrected.WE
		if(plot.raster==TRUE) {
			plot(corrected.WE_raster, main="Corrected Weighted Endemism")
		} #cls if(plot...
		outputs <- list(CWE = corrected.WE, corrected.WE_raster = corrected.WE_raster, weights = ranges, grid.matrix = cell_occur_matrix)
		return(outputs)
	} #cls if(type = c...


} #cls function

