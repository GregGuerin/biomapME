phylogenetic.endemism <- function(species_records, records="single", site.coords, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", sep.comm.spp = " ", phylo.tree, sep.phylo.spp = "_", frame.raster, deg.resolution=c(0.25,0.25), extent.vector, pe.type="weighted", plot.raster=TRUE, weight.type="cell", outlier_pct=95, own.weights, own.grid.matrix, own.phylo.cell.matrix, own.phyloMatrix, pd.standard=FALSE)
{


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

	sp::coordinates(species_records) <- c("LONGITUDE", "LATITUDE")

	if(class(phylo.tree) != "phylo") {
		stop("Phylogenetic tree must be in 'phylo' format")
	}

	if(missing(frame.raster)) {
		if(!(missing(own.grid.matrix))) {stop("You must supply a frame.raster with cells that match own.grid.matrix")}
		if(!(missing(own.phylo.cell.matrix))) {stop("You must supply a frame.raster with cells that match own.phylo.cell.matrix")}
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


	if(pe.type == "unweighted") {


		if(!(missing(own.grid.matrix))) {
		cat("Reading user-defined gridded occurrence matrix", "\n")
		if(class(own.grid.matrix) != "data.frame") {stop("Supplied own.grid.matrix must be a data.frame")}
		cell_occur_matrix <- own.grid.matrix
		} #clse if(!(missing(own.grid...



		if(missing(own.grid.matrix)) {

			if(sep.comm.spp != "none") {
				cat("Re-formatting species names", "\n") #this is to ensure mama() does not replace spaces with dots between genus/species
				removethis <- strsplit(as.character(species_records$SPECIES), sep.comm.spp, fixed=TRUE)
				removethis <- unlist(lapply(removethis, function(x) paste(x[1], x[2], sep="_", fixed=TRUE)))
				species_records$SPECIES <- factor(removethis)
			} #close if(sep.comm.spp...

			cat("Generating the gridded occurrence matrix", "\n")
			cell_numbers <- cellFromXY(frame.raster, species_records)
			cell_occur_matrix_prep <- data.frame(cell=cell_numbers, species=species_records$SPECIES, presence=rep(1, length(cell_numbers)))
			cell_occur_matrix_prep$species <- factor(cell_occur_matrix_prep$species)
			if(any(duplicated(cell_occur_matrix_prep))) {
				cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(duplicated(cell_occur_matrix_prep)),]
			} #cls if(any(duplicated))
			if(any(is.na(cell_occur_matrix_prep$cell))) {
				cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(is.na(cell_occur_matrix_prep$cell)),]
			} #cls if(any(is.na))
			cell_occur_matrix <- simba::mama(cell_occur_matrix_prep)
			}#cls if own matrix missing


		if(sep.comm.spp != "none") {
			removethis <- strsplit(as.character(colnames(cell_occur_matrix)), "_", fixed=TRUE)
			removethis <- unlist(lapply(removethis, function(x) paste(x[1], x[2], sep=sep.phylo.spp)))
			colnames(cell_occur_matrix) <- factor(removethis)
		} #close if(sep.comm.spp...

		cat("Occurrence matrix generated with dimensions: ", dim(cell_occur_matrix), "\n")



		cat("Checking that occurrence matrix and phylo match." , "\n")
		#trim excess species in cell_occur_matrix that aren't in phylo.tree
		cell_occur_matrix <- cell_occur_matrix[,which(colnames(cell_occur_matrix) %in% phylo.tree$tip.label)]
		#trim excess tips in phylo.tree that aren't in cell_occur_matrix
		phylo.tree <- drop.tip(phylo.tree, which(!(phylo.tree$tip.label %in% colnames(cell_occur_matrix))))



		cat("Calculating phylogenetic diversity", "\n")

		if(!pd.standard) {
		  
		PDmatFunc <- function(x, phylo.tree) {
			x <- x[,which(x>0)]
			drops <- setdiff(phylo.tree$tip.label, colnames(x))
			if(class(x) == "data.frame") {sub <- drop.tip(phylo.tree, drops)
			return(sum(sub$edge.length))} else {return(NA)}
			}

		PD <- list()
		for(i in 1:nrow(cell_occur_matrix)) {
			cell_row <- cell_occur_matrix[i,]
			PD[i] <- PDmatFunc(x=cell_row, phylo.tree=phylo.tree)
			}
		PD <- unlist(PD)
		names(PD) <- row.names(cell_occur_matrix)
		PD <- PD[!is.na(PD)]

		} #end if pd.standard=F

		
		if(pd.standard) {
		  cat("Standardising phylogenetic diversity by species richness", "\n")
		  PD <- PhyloMeasures::pd.query(tree=phylo.tree, matrix=cell_occur_matrix, standardize=TRUE)
		  names(PD) <- row.names(cell_occur_matrix)
		  PD <- PD[!is.na(PD)]
		} #end if pd.standard=T
		
		phyloDraster <- frame.raster
		phyloDraster[] <-  NA
		phyloDraster[as.numeric(as.character(names(PD)))] <- PD


		if(plot.raster == TRUE) {
			dev.new()
			plot(phyloDraster, main = "Phylogenetic diversity")
		} #cls if(plot.raster...


		cat("Collating outputs", "\n")
		output <- list()
		output$PD_raster <- phyloDraster
		output$PD <- PD
		output$grid.matrix <- cell_occur_matrix
		return(output)

	} #cls if pe.type= unweighted



	if(pe.type == "weighted") {

		if(!(missing(own.phylo.cell.matrix))) {
			cat("Reading user-defined gridded phylogenetic occurrence matrix.", "\n")
			if(class(own.phylo.cell.matrix) != "matrix") {stop("Supplied own.phylo.cell.matrix must be a matrix.")}
			phylo_occur_matrix <- own.phylo.cell.matrix

			if(missing(own.phyloMatrix)) {stop("You must also supply own.phyloMatrix", "\n")}
			if(class(own.phyloMatrix) != "matrix") {stop("Supplied own.phyloMatrix must be a matrix")}
			cat("Reading user-defined phyloMatrix", "\n")
			phyloMatrix <- own.phyloMatrix

			if(missing(own.grid.matrix)) {stop("You must also supply own.grid.matrix.", "\n")
			} #clse if(!(missing(own.grid..
			if(class(own.grid.matrix) != "data.frame") {stop("Supplied own.grid.matrix must be a data.frame.")}
		cat("Reading user-defined gridded occurrence matrix", "\n")
		cell_occur_matrix <- own.grid.matrix

		} #close if(!(missing(own.phylo.cell.matrix)))



		if(missing(own.phylo.cell.matrix)) {

		if(!(missing(own.grid.matrix))) {
		cat("Reading user-defined gridded occurrence matrix", "\n")
		if(class(own.grid.matrix) != "data.frame") {stop("Supplied own.grid.matrix must be a data.frame")}
		cell_occur_matrix <- own.grid.matrix

		} #clse if(!(missing(own.grid...



		if(missing(own.grid.matrix)) {

			if(sep.comm.spp != "none") {
				cat("Re-formatting species names", "\n") #this is to ensure mama() does not replace spaces with dots between genus/species
				removethis <- strsplit(as.character(species_records$SPECIES), sep.comm.spp, fixed=TRUE)
				removethis <- unlist(lapply(removethis, function(x) paste(x[1], x[2], sep="_", fixed=TRUE)))
				species_records$SPECIES <- factor(removethis)
			} #close if(sep.comm.spp...

			cat("Generating the gridded occurrence matrix", "\n")
			cell_numbers <- cellFromXY(frame.raster, species_records)
			cell_occur_matrix_prep <- data.frame(cell=cell_numbers, species=species_records$SPECIES, presence=rep(1, length(cell_numbers)))
			cell_occur_matrix_prep$species <- factor(cell_occur_matrix_prep$species)
			if(any(duplicated(cell_occur_matrix_prep))) {
				cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(duplicated(cell_occur_matrix_prep)),]
			} #cls if(any(duplicated))
			if(any(is.na(cell_occur_matrix_prep$cell))) {
				cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(is.na(cell_occur_matrix_prep$cell)),]
			} #cls if(any(is.na))
			cell_occur_matrix <- simba::mama(cell_occur_matrix_prep)
			}#cls if own matrix missing


		if(sep.comm.spp != "none") {
			#match species name formatting to phylo.tree ("_"). NOTE again, assuming you want two words in species names only. Because mama() turns spaces into dots, and either comm or phy data could have spaces, must repeat this step.
			removethis <- strsplit(as.character(colnames(cell_occur_matrix)), "_", fixed=TRUE)
			removethis <- unlist(lapply(removethis, function(x) paste(x[1], x[2], sep=sep.phylo.spp)))
			colnames(cell_occur_matrix) <- factor(removethis)
		} #close if(sep.comm.spp...


		cat("Occurrence matrix generated with dimensions: ", dim(cell_occur_matrix), "\n")





		cat("Checking that occurrence matrix and phylo match." , "\n")
		#trim excess species in cell_occur_matrix that aren't in phylo.tree
		cell_occur_matrix <- cell_occur_matrix[,which(colnames(cell_occur_matrix) %in% phylo.tree$tip.label)]
		#trim excess tips in phylo.tree that aren't in cell_occur_matrix
		phylo.tree <- drop.tip(phylo.tree, which(!(phylo.tree$tip.label %in% colnames(cell_occur_matrix))))



		###Conversion to matrix representation steps have been modified from the 'phylo.endemism.R' function (GNU GPL: http://davidnipperess.blogspot.com.au/2012/07/phyloendemism-r-function-for.html) of David Nipperess.
		cat("Converting phylo.tree to matrix", "\n")
		edgeNames <- as.character(1:length(phylo.tree$edge.length))
		tip.names <- phylo.tree$tip.label
		phyloMatrix <- matrix(0, length(phylo.tree$tip.label), length(phylo.tree$edge.length), dimnames = list(a=tip.names, b=edgeNames))
		for (i in 1:length(phylo.tree$tip.label)) {
			lineage <- which(phylo.tree$edge[,2] == i)
			node <- phylo.tree$edge[lineage,1]
			while (node > length(phylo.tree$tip.label)+1) {
				branch <- which (phylo.tree$edge[,2] == node)
				lineage <- c(lineage, branch)
				node <- phylo.tree$edge[branch,1]
			} #cls while(node)
		phyloMatrix[i,lineage] = 1
		} #cls for (i in...


		cat("Re-ordering the matrices to match...", "\n")
		cell_occur_matrix <- cell_occur_matrix[ ,sort.list(colnames(cell_occur_matrix))]
		cell_occur_matrix <- as.matrix(cell_occur_matrix)
		phyloMatrix <- phyloMatrix[sort(row.names(phyloMatrix)), ]
		cat("Phylo matrix generated with dimensions, ", dim(phyloMatrix), "and range ", range(phyloMatrix), "\n")


		cat("Generating a community phylogenetic branch matrix.", "\n")
		phylo_occur_matrix <- cell_occur_matrix %*% phyloMatrix #matrix product to identify which branches occur in which map cell by species shared by branches and cells
		phylo_occur_matrix[phylo_occur_matrix > 0] <- 1 #convert to binary occurrences of branches in cells


		}#cls if(missing(own.phylo.cell.matrix))


		cat("Calculating range weights", "\n")

		if(!(missing("own.weights"))) {
			if(class(own.weights) != "numeric") {stop("Supplied branch weights must be numeric")}
			if(length(own.weights) != length(phylo.tree$edge.length)) {stop("Supplied branch weights vector is different length than number of branches in supplied tree - must supply a weighting for each branch.")}
			cat("Invoking user-supplied weights", "\n")

			inv_rang_phylo_occur_matrix <- phylo_occur_matrix
			for (i in 1:ncol(inv_rang_phylo_occur_matrix)) {inv_rang_phylo_occur_matrix[,i] <- inv_rang_phylo_occur_matrix[,i]/own.weights[i]}
			ranges <- own.weights
		} #cls if(!(missing...



		if(missing(own.weights)) {


		if(weight.type == "cell") {
			ranges <- colSums(phylo_occur_matrix) #simply the number of cells (rows) each edge (column) occurs in.
			inv_rang_phylo_occur_matrix <- apply(phylo_occur_matrix, 2, function(x) {x/sum(x)})
		}	#cls if(weight.type == "cell)




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

			cell_centroids <- as.data.frame(coordinates(frame.raster))
			colnames(cell_centroids) <- c("LONGITUDE", "LATITUDE")
			cell_centroids$cells <- row.names(cell_centroids)

			edge_ranges <- function(x) {
				v <- rep(0, ncol(x))
				for (i in 1:ncol(x)) {
					temp <- as.data.frame(x[,i])
					colnames(temp) <- "edge_i"
					row.names(temp) <- row.names(x)
					temp$delete <- temp$edge_i
					if(any(temp$edge_i == 0)) {
						temp <- temp[-which(temp$edge_i == 0),]
					} #cls if(any(temp...
					temp$cells <- row.names(temp)
					temp <- merge(temp, cell_centroids, by="cells")[,c("LONGITUDE", "LATITUDE")]
					sp::coordinates(temp) <- c("LONGITUDE", "LATITUDE")
					
					cell_dimensions <- (mean(values(raster::area(frame.raster)))*1000000)^(0.5) #area in km^2 so converting to m^2 to match areaPolygon, then find square root to get back to 1-dimension
					
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
								edge_i_range_polygon <- try(adehabitatHR::mcp(temp, percent=outlier_pct))
								if(class(edge_i_range_polygon)[1] == "try-error") {
									v[i] <- max(CalcDists(temp))
								} #cls if(class(edge...
								if(!class(edge_i_range_polygon)[1] == "try-error") {
									v[i] <- max(CalcDists(as.data.frame(edge_i_range_polygon@polygons[[1]]@Polygons[[1]]@coords)))
								}	 #cls if(!class(spp...
								if(v[i] < cell_dimensions) {
									v[i] <- cell_dimensions
									warning("Range is less than spatial grain of frame.raster, returning span of single grid cell for ", colnames(x[i]))
								} #cls if(v[i] <...
							}	#cls if(nrow(temp) > 4)...
						} #cls else...


					cat("Edge complete:", i, v[i], "\n")

				} #cls for (i in...

				return(v)
			} #cls edge_ranges function

			ranges <- edge_ranges(phylo_occur_matrix)

			inv_rang_phylo_occur_matrix <- phylo_occur_matrix
			for (i in 1:ncol(inv_rang_phylo_occur_matrix)) {inv_rang_phylo_occur_matrix[,i] <- inv_rang_phylo_occur_matrix[,i]/ranges[i]}
		} #close if(weight.type = geo...

		} #cls if(missing(own.weights))





	#Now general procedures for all weight types
	cat("Calculating phylogenetic endemism.", "\n")
	for (i in 1:nrow(inv_rang_phylo_occur_matrix)) {
		inv_rang_phylo_occur_matrix[i,] <- inv_rang_phylo_occur_matrix[i,]*phylo.tree$edge.length
	}

	pd <- rowSums(inv_rang_phylo_occur_matrix)

	pd_raster <- frame.raster
	pd_raster[] <-  NA
	pd_raster[as.numeric(as.character(names(pd)))] <- pd
	if(plot.raster == TRUE) {
		dev.new()
		plot(pd_raster, main = "Phylogenetic endemism")
	} #cls if(plot...
	names(ranges) <- colnames(phyloMatrix)
	edge.lengths <- phylo.tree$edge.length
	names(edge.lengths) <- colnames(phyloMatrix)
	outputs <- list(PD = pd, ranges = ranges, PD_raster = pd_raster, phylo.cell.matrix = phylo_occur_matrix, edge.lengths = edge.lengths, phyloMatrix = phyloMatrix, grid.matrix = as.data.frame(cell_occur_matrix))

	return(outputs)
		} #cls if(pe.type == "weighted") {
	} #cls function
