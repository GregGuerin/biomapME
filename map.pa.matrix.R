map.pa.matrix <- function(species_records, records="single", site.coords, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", frame.raster, deg.resolution=c(0.25,0.25), extent.vector)
#
#Description --
#
#Given georeferenced incidence data for species, generates a binary presence/absence matrix associated with grid cells of a raster.
#
#Usage --
#
#For example:
#####Preparation for this example:#require(vegan) #load the vegan package#data(mite)#data(mite.xy) #load datasets from vegan####Usage of the function:#map.pa.matrix(mite, records="site", site.coords=mite.xy)
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
#frame.raster: an existing raster object. User can elect to supply a raster, in which case presences and absences are scored for grid cells in the raster. If not specified, the function generates a raster based on default or specified extent and resolution.
#
#deg.resolution/extent.vector: arguments specifying the map resolution and extent in degrees the user wishes the calculations and mapping to use. If none are specified, default resolution (0.25) and extent (data extent) are used. If a frame.raster is specified, these arguments are ignored.
#
#Details --
#
#This function generates a binary species presence/absence matrix associated with a raster layer based on georeferenced incidence data. This is a data processing step for mapping various biodiversity metrics onto raster layers. The outputs can be used as inputs into these functions, or if desired they can be used like site-based data (at the resolution of the raster) for various analysis such as ordination, or incidence/frequency data for particular species can be extracted.
#
#Value --
#
#Returns a list of length 2:
#
#$grid.matrix : a binary data.frame of species occurrences against grid cell numbers that match those in $pa.raster.
#
#$pa.raster : a ‘raster’ object for which species presence/absence is scored in $grid.matrix.
#
#Required packages --
#
#simba, raster
#
#Authors --
#
#Greg R. Guerin
#
#References --
#
#Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12361
#
#License --
#
#GNU GPL-3
#
{
	require(raster)
	require(simba)
	
		
	if(class(species_records) != "data.frame") {
		stop("Species data must be in a data.frame")
	} #cls if(class(species…
	
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
	} #cls if(records == "single. . .
	
	
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
	} #cls if(missing. . .
	
	
	if(!(extent(species_records)@xmin >= extent(frame.raster)@xmin & extent(species_records)@xmax <= extent(frame.raster)@xmax & extent(species_records)@ymin >= extent(frame.raster)@ymin & extent(species_records)@ymax <= extent(frame.raster)@ymax)) {
		cat("Some point locations lie outside the frame raster -- trimming these records", "\n")
		species_record_COORDS <- as.data.frame(coordinates(species_records))
		species_records <- species_records[-which(species_record_COORDS$LONGITUDE < extent(frame.raster)@xmin | species_record_COORDS$LONGITUDE > extent(frame.raster)@xmax | species_record_COORDS$LATITUDE < extent(frame.raster)@ymin | species_record_COORDS$LATITUDE > extent(frame.raster)@ymax),]
		} #cls if(!(extent(. . .
				

			
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

		
		
				
		cat("Occurrence matrix generated with dimensions: ", dim(cell_occur_matrix), "\n")
		
		
		pa.raster <- frame.raster
		pa.raster[] <-  NA
			
				
		cat("Collating outputs", "\n")
		output <- list()
		output$pa.raster <- pa.raster
		output$grid.matrix <- cell_occur_matrix
		return(output)
				

	} #cls function