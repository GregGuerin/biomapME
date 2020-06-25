convert.site.data <- function(species_records, site.coords)
{
  if(class(species_records) != "data.frame") {
    stop("Species data must be in a data.frame")
  }
  if(class(site.coords) != "data.frame") {
    stop("Spatial coordinates must be in a data.frame")
  }

  dat <-  data.frame(SPECIES = "hold", LONGITUDE = 0, LATITUDE = 0)
  nam <-  names(species_records)
  for(ii in 1:ncol(species_records)){
    w <-  species_records[,ii]>0
    dat <-  rbind(dat, setNames(data.frame(rep(nam[ii],sum(w)),site.coords[w,]), names(dat)))
  }

  return(droplevels(dat[-1,]))

}




