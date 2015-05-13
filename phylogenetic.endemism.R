phylogenetic.endemism <- function(species_records, records="single", site.coords, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", sep.comm.spp = " ", phylo.tree, sep.phylo.spp = "_", frame.raster, deg.resolution=c(0.25,0.25), extent.vector, pe.type="weighted", plot.raster=TRUE, weight.type="cell", outlier_pct=95, own.weights, own.grid.matrix, own.phylo.cell.matrix, own.phyloMatrix)

#shorter version, showing relevant arguments for inventory style data (individual records):
#phylogenetic.endemism <- function(species_records, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", sep.comm.spp = " ", phylo.tree, sep.phylo.spp = "_", ...)
#
#shorter version, showing relevant arguments for site/plot based data:
#phylogenetic.endemism <- function(species_records, records="single", site.coords, sep.comm.spp = " ", phylo.tree, sep.phylo.spp = "_", ...)
#
#Description --
#
#Calculates phylogenetic endemism (phylogenetic diversity inversely weighted by the spatial range of particular branch lengths) or alternatively (unweighted) phylogenetic diversity across gridded maps using individual or site-based point records.
#
#Usage --
#
#For example:
#####Preparation for this example:#library(vegan) #load the vegan package#data(mite)#data(mite.xy) #load datasets from vegan#library(ape) #for rtree function:#mite.tree <- rtree(n=ncol(mite), tip.label=colnames(mite)) #for this example, generate a phylogenetic tree of the species in the mite dataset with random relationships and branch lengths####Usage of the function:
#source(“phylogenetic.endemism.R”)#mite.PE <- phylogenetic.endemism(mite, records="site", site.coords=mite.xy, sep.comm.spp="none", phylo.tree=mite.tree, sep.phylo.spp="none", weight.type="geo")
#
#Arguments --
#
#species_records: a data.frame, either with:
#a) rows as individual species records, and columns that include fields for species name, longitude and latitude (see ‘species’, ‘longitude’, ‘latitude’ below); or
#b) rows as sites and columns as species, in which case ‘site.coords’ (below) must also be supplied
#
#records: are the species_records in single/long format (the default, records="single") or in site-based/short format (records="site")
#
#site.coords: for site-based data (records="site"), a data.frame with rows as sites (/field plots) (names match the row names of species_records) and their geographic longlat coordinates: first column must be x/longitude, second column y/latitude
#
#species: for records="single" (i.e. individual occurrence data); what colname in the supplied species_records contains species names?
#
#latitude: for records="single"; what colname in the supplied species_records contains latitude values?
#
#longitude: for records="single"; what colname in the supplied species_records contains longitude values?
#
#sep.comm.spp: the genus_species separating character in the community data. If there is none (i.e. taxon names represented by a single 'word'), argument should be set to sep.comm.spp=“none". The purpose of this argument is that it is common for separators to differ between community data and tree due to the different processing functions applied, and to avoid having to reformat the data, this argument allows them to be matched within the function. Default is a space.
#
#phylo.tree: a phylogenetic tree of class 'phylo' containing species in the occurrence data. Must have branch lengths.
#
#sep.phylo.spp: the genus_species separating character in the phylo.tree (tip.labels). If there is none (i.e. taxon names represented by a single 'word'), argument should be set to sep.phylo.spp=“none". See sep.comm.spp argument for purpose. Default for phylo.tree is an underscore, as this is a common format in community ecology using e.g. phylomatic trees.
#
#frame.raster: an existing raster object. User can elect to supply a raster as the frame for calculations and mapping. If not specified, the function will generate a raster based on default or specified extent and resolution.
#
#deg.resolution/extent.vector: arguments specifying the map resolution and extent in degrees the user wishes the calculations and mapping to use. If none are specified, default resolution (0.25) and extent (data extent) are used. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster).
#
#pe.type: either "weighted" (default; phylogenetic endemism), or "unweighted" (pyhlogenetic diversity).
#
#plot.raster: whether or not to plot the output raster of endemism scores. Whether plotted or not, a raster object is stored in the output.
#
#weight.type: default is "cell" (cell-based range weights), while "geo" will calculate georeferenced range 'span' weights. Argument is ignored if pe.type="unweighted" or if own.weights are provided.
#
#outlier_pct: for the calculation of range 'span' via convex polygons, this argument can be used to remove outliers on a percentage basis via the 'mcp' function in package adehabitat, i.e. 95 (default) means 5% most outlying points are removed.
# 
#own.weights: an optional user-supplied numeric vector of branch length weights for calculating phylogenetic endemism. This argument is intended mainly so that time consuming calculation of georeferenced weights can be done once and the result stored and used for subsequent re-runs. Alternatively, user-calculated weights can be used in subsequent runs. The weights must be matched to a matrix representation of phylo.tree (see below).
#
#own.grid.matrix: user supplied binary matrix of species against grid cell numbers, rather than this being generated within the function. The purpose of this argument is to save time for repeat runs, given the step can be time consuming for large datasets. Must be provided if own.phylo.cell.matrix is supplied-- although the function no longer requires the grid.matrix, it is still needed for inclusion in the outputs, so that they can be subsequently fed into the pe.null.test.R function (which does requires the grid.matrix...). If supplied, a frame.raster must also be supplied that has cell numbers which match the row.names of the own.grid.matrix. Assumes that genus/species words are separated by the sep.phylo.spp, as this is how the grid.matrix is returned from the function originally
#
#own.phylo.cell.matrix : user can supply the matrix representing occurrences of phylogenetic branches in map grid cells from previous runs to save computation time (i.e. $phylo.cell.matrix of a previous phylogenetic.endemism.R output). A matching 'frame.raster' must also be supplied.
#
#own.phyloMatrix : if user provides the above own.phylo.cell.matrix in a repeat run of the function, this must also be supplied, as the function still defaults to returning 'phylo.matrix' even though the function now skips its generation. e.g. $phyloMatrix of an previous phylogenetic.endemism.R run).
#
#Details --
#
#This implementation of phylogenetic endemism allows alternative calculation of weights for branch length ranges. Weights can be calculated based on the frequency of occurrence in grid cells, or alternatively by the georeferenced span of the range. Unweighted phylogenetic diversity can also be selected.
#
#Value --
#
#Returns a list of length 3 (pe.type="unweighted") or 7 (pe.type="weighted"):
#
#PE and PD:
#$PD : vector of phylogenetic diversity/endemism scores.
#
#$PD_raster: raster map with phylogenetic diversity/endemism scores.
#
#$grid.matrix : a binary data.frame of species against grid cell numbers used in the function which is returned so that it can be re-used to save computation time, and because it is required in the downstream pe.null.test.R function to calculated species sample probabilities. Note that this only includes species that are found in the phylo.tree
#
#PE only:
#$ranges : a named numeric vector of weights used to calculate endemism (equivalent to range size in metres if weight.type="geo", range size in cells if weight.type="cell" (default), or the user supplied weights if own.weights was supplied) (skipped if pe.type="unweighted").
#
#$phyloMatrix : a binary matrix representation of the phylo.tree, returned mainly for convenience for downstream use (skipped if pe.type="unweighted").
#
#$phylo.cell.matrix : a binary matrix recording the presence of particular phylogenetic branches in map grid cells used in the function, which is returned for re-use in subsequent runs for efficiency (e.g. with different weights), (skipped if pe.type="unweighted").
#
#$edge.lengths : a numeric vector of edge lengths from phylo.tree
#
#Required packages --
#
#simba, geosphere, adehabitat, raster, ape
#
#Authors --
#
#Greg R. Guerin
#
#References --
#
#Guerin, G.R. and Lowe, A.J. (submitted) Mapping phylogenetic endemism in R using georeferenced branch extents.
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
	require(raster)
	require(simba)
	require(adehabitat)
	require(geosphere)
	require(ape)
	
	
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
				dat <- rbind(dat, setNames(data.frame(rep(nam[ii],sum(w)),site.coords[w,]), names(dat)))
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
			cell_occur_matrix <- mama(cell_occur_matrix_prep)
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
			cell_occur_matrix <- mama(cell_occur_matrix_prep)
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
					temp <- merge(temp, cell_centroids, by="cells")
					temp <- data.frame(LONGITUDE=temp$LONGITUDE, LATITUDE=temp$LATITUDE)
					cell_dimensions <- (mean(values(raster::area(frame.raster)))*1000000)^(0.5) #area in km^2 so converting to m^2 to match areaPolygon, then find square root to get back to 1-dimension
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
								edge_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
								if(class(edge_i_range_polygon)[1] == "try-error") {
									v[i] <- max(CalcDists(temp))
								} #cls if(class(edge...
								if(!class(edge_i_range_polygon)[1] == "try-error") {
									v[i] <- max(CalcDists(as.data.frame(edge_i_range_polygon[,2:3])))
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