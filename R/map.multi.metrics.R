map.multi.metrics <- function(species_records, site.coords, alpha=TRUE, beta=TRUE, frame.raster, deg.resolution=c(0.25,0.25), extent.vector, plot.raster=TRUE) {
  
  if(class(species_records) != "data.frame") {
    stop("Species data must be in a data.frame")
  }
  
  if(missing(frame.raster)) {
    frame.raster <- raster::raster()
    frame.raster[] <- NA
    if(missing(extent.vector)) {
      raster::extent(frame.raster)@xmin <- floor(min(site.coords[,1]))
      raster::extent(frame.raster)@ymin <- floor(min(site.coords[,2]))
      raster::extent(frame.raster)@xmax <- ceiling(max(site.coords[,1]))
      raster::extent(frame.raster)@ymax <- ceiling(max(site.coords[,2]))
    }
    if(!(missing(extent.vector))) {
      extent(frame.raster) <- extent.vector
    }
    res(frame.raster) <- deg.resolution
    cat("Generating frame raster at ", deg.resolution, " resolution and extent defined by: ", extent(frame.raster)@xmin, extent(frame.raster)@xmax, extent(frame.raster)@ymin, extent(frame.raster)@ymax,"\n")
  }
  
  site.coords$cell <- raster::cellFromXY(frame.raster, site.coords)
  
###########################################
  if(beta) {
  
    multi_list <- split(species_records, site.coords$cell) #list of occurrence matrices split by shared grid cells
    
###
  beta.multi.cell <- function(comm) {
    try(BAT::beta.multi(comm))
  }
  ###
  beta_result <- lapply(multi_list, beta.multi.cell) #get beta metrics for each matrix
  
  

  

  total_mean <- lapply(beta_result, function(x) return(x[1])) #extract first value from each table in the list (it is a list too) - equals total beta mean
  
  success_cells <- names(total_mean[unlist(lapply(total_mean, is.numeric))]) #cell numbers where there was a result (data in cell and metric calculated, i.e. >1 plot)
  
  beta_result <- as.data.frame(do.call(rbind, lapply(beta_result, function(x) {return(t(x)[1,])})[success_cells]))
  

  
frame.raster[as.numeric(rownames(beta_result))] <- beta_result$Btotal #add worked values to matching cells
  
 results <- list()
  results$beta_result <- beta_result
  results$beta_raster <- frame.raster
  
  } #end if beta  
  
  #####################################################
  
  if(alpha) {
    
    alpha_result <- vegan::specpool(species_records, site.coords$cell)
    
    frame.raster[] <- NA
    frame.raster[as.numeric(rownames(alpha_result))] <- alpha_result$chao
    
    results$alpha_result <- alpha_result
    results$alpha_raster <- frame.raster
  } #end if alpha

  
##############################################
 return(results)
  
}