map.multi.metrics <- function(species_records, site.coords, alpha=TRUE, beta=TRUE, phylogenetic=FALSE, phylo.tree, frame.raster, deg.resolution=c(0.25,0.25), extent.vector, custom_grouping, plot.raster=TRUE) {
  
  
  #checks
  if(class(species_records) != "data.frame") {
    stop("Species data must be in a data.frame")
  }
  
  if(class(site.coords) != "data.frame") {
    stop("Species data must be in a data.frame")
  }
  
  if(phylogenetic) {
    if(missing(phylo.tree)) {
      stop("You must provide a phylo.tree if phylogenetic = TRUE")
    }
    if(class(phylo.tree) != "phylo") {
      stop("Phylogenetic tree must be in 'phylo' format")
    }
  }
  
  if(!missing(custom_grouping)) {
    if(class(custom_grouping) != "numeric") {
      stop("Argument custom_grouping must be a numeric constant indicating a column of site.coords")
    }
    if(!is.vector(site.coords[ , custom_grouping])) {
      stop("Cannot locate grouping column in site.coords")
    }
    cat("Calculating metrics by supplied site groupings - no raster will be produced.", "\n")
    site.coords$cell <- site.coords[ , custom_grouping]
  }
  
  if(any(!is.logical(c(alpha, beta, phylogenetic, plot.raster)))) {
    stop("Arguments alpha, beta, phylogenetic and plot.raster must be TRUE or FALSE")
  }
  
  if(!alpha & phylogenetic) {
    stop("Phylogenetic diversity estimation relies on alpha diversity. Set alpha = TRUE to calculate PD")
  }
  
  if(!missing(frame.raster)) {
    if(class(frame.raster != "RasterLayer")) {
      stop("frame.raster must be a RasterLayer object")
    }
  }
  
  if(missing(custom_grouping)) {
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
  }
  
  
########################

results <- list()

###########################################

  if(alpha) {
    
    alpha_result <- vegan::specpool(species_records, site.coords$cell)
    results$alpha_result <- alpha_result
    
    if(missing(custom_grouping)) {
      frame.raster[] <- NA
    frame.raster[as.numeric(rownames(alpha_result))] <- alpha_result$chao
    results$alpha_raster <- frame.raster
    }
    
  } #end if alpha
  
  ##############################################
  

  if(beta) {
  
    multi_list <- split(species_records, site.coords$cell) #list of occurrence matrices split by shared grid cells
    
###
  beta.multi.cell <- function(comm) {
    try(BAT::beta.multi(comm), silent=TRUE)
  }
  ###
  beta_result <- lapply(multi_list, beta.multi.cell) #get beta metrics for each matrix
  
  total_mean <- lapply(beta_result, function(x) return(x[1])) #extract first value from each table in the list (it is a list too) - equals total beta mean
  
  success_cells <- names(total_mean[unlist(lapply(total_mean, is.numeric))]) #cell numbers where there was a result (data in cell and metric calculated, i.e. >1 plot)
  
  beta_result <- as.data.frame(do.call(rbind, lapply(beta_result, function(x) {return(t(x)[1,])})[success_cells]))
  
  results$beta_result <- beta_result
  

  if(missing(custom_grouping)) {
    frame.raster[] <- NA
    frame.raster[as.numeric(rownames(beta_result))] <- beta_result$Btotal #add worked values to matching cells
    results$beta_raster <- frame.raster
  }
  
  } #end if beta  
  
  #####################################################
  
  
  if(phylogenetic) {
    
    cat("Checking that occurrence matrix and phylo match." , "\n")
    #trim excess species in species_records that aren't in phylo.tree
    species_records <- species_records[,which(colnames(species_records) %in% phylo.tree$tip.label)]
    #trim excess tips in phylo.tree that aren't in species_records
    phylo.tree <- ape::drop.tip(phylo.tree, which(!(phylo.tree$tip.label %in% colnames(species_records))))
    
    if(!exists("multi_list")) {
      multi_list <- split(species_records, site.coords$cell) #list of occurrence matrices split by shared grid cells
    }
    
      cat("Calculating phylogenetic diversity", "\n")
    
    PDmatFunc <- function(x, phylo.treeee=phylo.tree) {
      x <- x[x > 0]
      sub <- ape::drop.tip(phylo.treeee, which(!(phylo.treeee$tip.label %in% names(x))))
      return(sum(sub$edge.length))
    }
    
    
    PDobs <- unlist(lapply(lapply(multi_list, colSums), PDmatFunc))
    c <- alpha_result$Species / alpha_result$chao #vector of sample completeness index
    PDest <- PDobs / c
    
    phylogenetic_result <- data.frame(PDobs=PDobs, completeness=c, PDest=PDest)
    rownames(phylogenetic_result) <- names(PDobs)
    results$phylogenetic_result <- phylogenetic_result
    
    if(missing(custom_grouping)) {
      frame.raster[] <- NA
      frame.raster[as.numeric(rownames(phylogenetic_result))] <- phylogenetic_result$PDest
      results$phylogenetic_raster <- frame.raster
    }
    
    
  } #end if phylogenetic
  
#############################################

if(missing(custom_grouping)) {
  if(plot.raster) {
    try(lapply(results[grep("raster", names(results))], plot))
  }
}

 return(results)
  
}