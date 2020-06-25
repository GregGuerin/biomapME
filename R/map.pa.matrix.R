map.pa.matrix <- function(species_records, records="single", site.coords, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", frame.raster, deg.resolution=c(0.25,0.25), extent.vector)
{

	if(class(species_records) != "data.frame") {
		stop("Species data must be in a data.frame")
	} #cls if(class(species…

	if(records == "site") {
		
		species_records <- convert.site.data(species_records, site.coords)
		
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
