weighted.endemism <- function(species_records, records="single", site.coords, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", frame.raster, deg.resolution=c(0.25,0.25), extent.vector, type="weighted", plot.raster=TRUE, own.weights, weight.type="cell", geo.type="cell", geo.calc="max.dist", outlier_pct=95, verbose=TRUE, own.grid.matrix)
#
#Description --
#
#Calculates (taxonomic / species) weighted endemism (species richness inversely weighted by species ranges) across gridded maps using single or site-based point records.
#
#Usage --
#
#For example:
#
#require(vegan)
#
#endemism_mydata <- weighted.endemism(mite, site.coords = mite.xy, records="site")
#
#endemism_mydata2 <- weighted.endemism(mite, site.coords = mite.xy, records="site", weight.type="geo", own.grid.matrix = endemism_mydata$grid.matrix, frame.raster=endemism_mydata$WE_raster)
#
#endemism_mydata3 <- weighted.endemism(mite, site.coords = mite.xy, records="site", own.weights = endemism_mydata2$weights)
#
#Arguments --
#
#species_records: a data.frame, either with:
#a) rows as individual species records, and columns that include fields for species name, longitude and latitude (see species, longitude, latitude below); or
#b) rows as sites and columns as species, in which case site.coords (below) must also be supplied
#
#records: are the species_records in single/long format (the default, "single") or in site-based/short format (records="site")
#
#site.coords: for site-based data (records="site"), a data.frame with the sites (/field plots) that match the column names of species_records and their longlat coordinates
#
#species: for records="single"; what colname in the supplied species_records contains species names?
#
#latitude: for records="single"; what colname in the supplied species_records contains latitude values?
#
#longitude: for records="single"; what colname in the supplied species_records contains longitude values?
#
#frame.raster: an existing raster object the user can optionally elect to supply as the frame for calculations and mapping
#
#deg.resolution/extent.vector: arguments specifying the map resolution and extent in degrees the user wishes the calculations and mapping to use. If no frame is specified, an arbitrary resolution and extent is supplied. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster)
#
#type: either "weighted" (default; weighted endemism), or "corrected" (corrected weighted endemism - the 'per-species' weighted endemism as per Crisp et al. (2001). This is provided for convenience, but is not particularly recommended - insterad, use the outputs of weighted.endemism in endemism.null.test function to compare to null expectations of endemism given the species richness)
#
#plot.raster: whether or not to plot the output raster with endemism scores. Either way, the raster object is stored in the output.
#
#own.weights: an optional user-supplied numeric vector of species weights for calculating endemism. Values must have names that are a complete and exact match for the species names in species_records as each species must have a weight. This optional argument is intended mainly so that the more time consuming calculation of geographic weights can be done once and the result stored and used for subsequent re-runs of the endemism calculations (see example under Usage above).
#
#weight.type: default is "cell" (cell-based range weights), while "geo" will calculate geographic range weights. Weight.type "richness" sets the weights to 1, which is equivalent to calculating simple species richness (note setting type="corrected" in this case will give each cell a score of 1). Argument is ignored if own.weights is supplied.
#
#geo.type: default is "cell". This argument is only considered when weight.type is set to "geo" (above), in which case geo.type="cell" calculates geographic ranges based on map grid cell centroids. This can optionally be set to geo.type="point", in which case the geographic range weightings are calculated based on the point locations of each species - this is supplied for reference, but is not especially recommended as it is slower and doesn’t provide much more information at the resolution and scale of most analyses. 
#
#geo.calc: default is "max.dist". This argument is only considered when when weight.type is set to "geo" and applies to both geo.type="cell" and geo.type="point". If geo.calc="max.dist" the default maximum geographic distance (‘span’) between cell centroids/point locations is calculated for a species range. If geo.calc is set to "polygon", the function weights species ranges based on the area of a MCP (minimum convex polygon) that contains all points. Further arguments to this function can be included, such as changing the default outlier_pct=95 (removes outlying locations). Additionally, the geo.type="point" option is not recommended for the "polygon" method, as it is more likely to lead to errors where nearby point locations do not allow drawing of a spanning polygon. In this case, cell-centroid based calculations ensure that multiple records are spatially separated (in different cells) and that occurrences within a single cell are returned as the area of that cell.
#
#outlier_pct: for the calculation of range span or area via convex polygons (at least 5 records of the species), this argument can be used to remove outliers on a percentage basis via the mcp function in package adehabitat, i.e. 95 means 5% most outlying points are removed. Default is 95.
#
#own.grid.matrix: user can supply a binary matrix of species against grid cell numbers, rather than this being generated within the function. The purpose of this argument is that the step can be time consuming for large datasets, so the user can return the matrix that is returned from the function in subsequent runs with different setting to improve speed. If this is supplied, a frame.raster must also be supplied that has cell numbers which match the row.names of the own.grid.matrix
#
#Details --
#
#This implementation of weighted endemism allows alternative calculation of weights for species ranges as well as the option of user-supplied weights. Weights can be calculated based on the frequency of occurrence in grid cells, or alternatively by the geographical size of the species range, calculated in one (span) or two (area) dimensions.
#
#Value --
#
#Returns a list of length 4:
#
#$WE (/$CWE) : vector of weighted endemism scores
#
#$WE_Raster (/$CWE) : raster map with endemism scores
#
#$weights : a named numeric vector of weights used to calculate endemism (equivalent to range size in metres if weight.type="geo", range size in cells if weight.type="cell" (default), or the user supplied weights if own.weights was supplied)
#
#$grid.matrix : a binary data.frame of species against grid cell numbers used in the function which is returned so that it can be re-used in subsequent runs to save time
#
#Required packages --
#
#simba, geosphere, adehabitat, raster
#
#Authors --
#
#Greg R. Guerin & Lasse Ruokolainen
#
#References --
#
#Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12361
#
#License --
#
#GPL-3
#
{
	
	require(raster)
	require(simba)
	require(adehabitat)
	require(geosphere)
	
	if(outlier_pct > 99 | outlier_pct < 1) {
		stop("Outlier_pct should be a percentage")
	}
	
	if(class(species_records) != "data.frame") {
		stop("Species data must be in a data.frame")
	}
	
	if(records == "site") {
		convert <- function(an.occurrence.matrix, site.coords) {
			dat <-  data.frame(SPECIES = "hold",LONGITUDE = 0,LATITUDE = 0)
			nam <-  names(an.occurrence.matrix)
			for(ii in 1:ncol(an.occurrence.matrix)){
				w <-  an.occurrence.matrix[,ii]>0
				dat <-  rbind(dat, setNames(data.frame(rep(nam[ii],sum(w)),site.coords[w,]), names(dat)))
			}
			return(dat[-1,])
		}
		species_records <- convert(species_records, site.coords)
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
		cat("Generating the gridded occurrence matrix", "\n")
		cell_numbers <- cellFromXY(frame.raster, species_records)
		cell_occur_matrix_prep <- data.frame(cell=cell_numbers, species=species_records$SPECIES, presence=rep(1, length(cell_numbers)))
		cell_occur_matrix_prep$species <- factor(cell_occur_matrix_prep$species)
		if(any(duplicated(cell_occur_matrix_prep))) {
			cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(duplicated(cell_occur_matrix_prep)),]
		}
		if(any(is.na(cell_occur_matrix_prep$cell))) {
			cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(is.na(cell_occur_matrix_prep$cell)),]
		}
		cell_occur_matrix <- mama(cell_occur_matrix_prep)
		cat("Occurrence matrix generated with dimensions: ", dim(cell_occur_matrix), "\n")
	} #cls if(missing(own.grid...
		
	
	
	
	if(!(missing("own.weights"))) {
		if(class(own.weights) != "numeric") {stop("Supplied species weights must be numeric")}
		if(!(all(names(own.weights) %in% colnames(cell_occur_matrix)) & all(colnames(cell_occur_matrix) %in% names(own.weights)))) {stop("Species names for supplied weights are not a complete match for species in supplied records")} #true if all in it both ways
		if(length(own.weights) != length(colnames(cell_occur_matrix))) {stop("Supplied weights vector is different length than number of species - must supply a weighting for each species")}
		
		cat("Calculating user supplied weights", "\n")
		inv_rang_cell_occur_mat <- cell_occur_matrix
		for (i in 1:ncol(inv_rang_cell_occur_mat)) {inv_rang_cell_occur_mat[,i] <- inv_rang_cell_occur_mat[,i]/own.weights[which(names(own.weights) == colnames(inv_rang_cell_occur_mat)[i])]}
		ranges <- own.weights
	} #cls if(!(missing...
	
	
	
	if(missing("own.weights")) {
		
		
		if(weight.type=="cell") {
			cat("Calculating cell-based range weights", "\n")
			inv_rang_cell_occur_mat <- apply(cell_occur_matrix, 2, function(x) {x/sum(x)})
			ranges <- apply(cell_occur_matrix, 2, function(x) {sum(x)}) #can we just do colSums?
		} #cls if(weight.type = cell...
		
		if(weight.type=="richness") {
			inv_rang_cell_occur_mat <- cell_occur_matrix
			ranges <- apply(cell_occur_matrix, 2, function(x) {x = 1})
		} #cls if(w.t = richness...
		
		if(weight.type=="geo") {	
			
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
							temp <- merge(temp, cell_centroids, by="cells")
							temp <- data.frame(LONGITUDE=temp$LONGITUDE, LATITUDE=temp$LATITUDE)
							
							if(geo.calc == "max.dist") {
								cell_dimensions <- (mean(values(area(frame.raster)))*1000000)^(0.5) #area in km^2 so converting to m^2 to match areaPolygon, then find square root to get back to 1-dimension
								if(nrow(temp) < 2) {
									v[i] <- cell_dimensions
									warning("Only one record, returning the span of single grid cell for ", colnames(x[i]))
									} else {
										if(nrow(temp) < 5) {
											v[i] <- max(CalcDists(temp))
											if(v[i] < cell_dimensions) {
												v[i] <- cell_dimensions
												warning("Range is less than spatial grain of frame.raster, returning span of single grid cell for ", colnames(x[i]))
											} #cls if(v[i] < cell_dimensions)...
										} #cls if nrow(temp) < 5...
										if(nrow(temp) > 4) {
											spp_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
											if(class(spp_i_range_polygon)[1] == "try-error") {
												v[i] <- max(CalcDists(temp))
											} #cls if(class(spp...
											if(!class(spp_i_range_polygon)[1] == "try-error") {
												v[i] <- max(CalcDists(as.data.frame(spp_i_range_polygon[,2:3])))
											}	 #cls if(!class(spp...
											if(v[i] < cell_dimensions) {
												v[i] <- cell_dimensions
												warning("Range is less than spatial grain of frame.raster, returning span of single grid cell for ", colnames(x[i]))
												} #cls if(v[i] <...
										}	#cls if(nrow(temp) > 4)...						
									} #cls else...
							} #cls if geo.calc = max.dist...
							
							if(geo.calc == "polygon") {
								cell_size <- mean(values(area(frame.raster)))*1000000 #area in km^2 so convert to m^2 to match areaPolygon
								if(nrow(temp) < 5) {
									polygon_area <- try(areaPolygon(temp))
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
								if(nrow(temp) > 4) {
									spp_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
									if(class(spp_i_range_polygon)[1] == "try-error") {
										v[i] <- cell_size
										warning("Unable to compute polygon, returning approximate area of single grid cell for ", colnames(x[i]))
									} #cls if(class(spp...
									if(class(spp_i_range_polygon)[1] != "try-error") {
										plot(spp_i_range_polygon, main=paste(polygon_area))
										polygon_area <- try(areaPolygon(as.data.frame(spp_i_range_polygon[,2:3])))
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
						colnames(temp) <- colnames(x) #?necessary
						temp <- data.frame(LONGITUDE=temp$LONGITUDE, LATITUDE=temp$LATITUDE)
								
								if(geo.calc == "max.dist") {
									cell_dimensions <- (mean(values(area(frame.raster)))*1000000)^(0.5) #area in km^2 so convert to m^2 to match areaPolygon, then find square root to get back to 1-dimension
									if(nrow(temp) < 2) {
										v[n] <- cell_dimensions
										warning("Only one record, returning the size of a single grid cell (1-dimension) for ", i)
									} else {
											if(nrow(temp) < 5) {
												v[n] <- max(CalcDists(temp))
												if(v[n] < cell_dimensions) {
													v[n] <- cell_dimensions
													warning("Range span less than spatial grain of frame.raster, returning span of single grid cell for ", i)
												} #cls if(v[n] < cell_dimensions...
											} #cls if(nrow(temp) < 5)...
											if(nrow(temp) > 4) {
												spp_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
												if(class(spp_i_range_polygon)[1] == "try-error") {
													v[n] <- max(CalcDists(temp))
												} #cls if(class...
												if(!class(spp_i_range_polygon)[1] == "try-error") {
													v[n] <- max(CalcDists(as.data.frame(spp_i_range_polygon[,2:3])))
												} #cls if(!(class...
												if(v[n] < cell_dimensions) {
													v[n] <- cell_dimensions
													warning("Range span less than spatial grain of frame.raster, returning span of single grid cell for ", i)
												} #cls if(v[n] <...
											} #cls if(nrow(temp) > 4)...
										} #cls else...
								} #cls if(geo.calc == "max.dist")... [three open now...]
								
								if(geo.calc == "polygon") {
									cell_size <- mean(values(area(frame.raster)))*1000000 #returns size in m^2 (convert from km^2)
									if(nrow(temp) < 5) {
										polygon_area <- try(areaPolygon(temp))
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
									if(nrow(temp) > 4) {
										spp_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
										if(class(spp_i_range_polygon)[1] == "try-error") {
											v[n] <- cell_size
											warning("Unable to compute a polygon, returning approximate area of single grid cell for ", i)
										} #cls if(class(spp...
										if(class(spp_i_range_polygon)[1] != "try-error") {
											plot(spp_i_range_polygon, main=paste(polygon_area))
											polygon_area <- try(areaPolygon(as.data.frame(spp_i_range_polygon[,2:3])))
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

