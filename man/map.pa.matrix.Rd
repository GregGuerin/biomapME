\name{map.pa.matrix}
\alias{map.pa.matrix}
\title{Map Presence-absence Matrix
}
\description{Given georeferenced incidence data for species, generates a binary presence/absence matrix associated with grid cells of a raster.
}
\usage{
map.pa.matrix(species_records, records = "single", site.coords, species = "SPECIES",
longitude = "LONGITUDE", latitude = "LATITUDE", frame.raster, deg.resolution = c(0.25, 0.25),
extent.vector)
}
\arguments{
  \item{species_records}{A data.frame, either with:
  a) rows as individual species records, and columns that include fields for species name, longitude and latitude (see 'species', 'longitude', 'latitude' below); or
  b) rows as sites and columns as species, in which case ‘site.coords’ (below) must also be supplied.
}
  \item{records}{Are the species_records in single/long format (the default, records="single") or in site-based/short format (records="site")?
}
  \item{site.coords}{For site-based data (records="site"), a data.frame with rows as sites (/field plots) (names match the row names of species_records) and their geographic coordinates: first column must be x/longitude, second column y/latitude.
}
  \item{species}{For records="single" (i.e. individual occurrence data); what colname in the supplied species_records contains species names?
}
  \item{longitude}{For records="single"; what colname in the supplied species_records contains longitude values?
}
  \item{latitude}{For records="single"; what colname in the supplied species_records contains latitude values?
}
  \item{frame.raster}{An existing rasterLayer object. User can elect to supply a raster, in which case presences and absences are scored for grid cells in the raster. If not specified, the function generates a raster based on default or specified extent and resolution.
}
  \item{deg.resolution}{Arguments specifying the map resolution (in degrees) the user wishes the calculations and mapping to use. If none are specified, default resolution (0.25) is used. If a frame.raster is specified, this argument is ignored.
}
  \item{extent.vector}{Argument specifying the map extent (in degrees) the user wishes the calculations and mapping to use. If none are specified, default extent (data extent) is used. If a frame.raster is specified, this argument is ignored.
}
}
\details{This function generates a binary species presence/absence matrix associated with a raster layer based on georeferenced incidence data. This is a data processing step for mapping various biodiversity metrics onto raster layers. The outputs can be used as inputs into these functions, or if desired they can be used like site-based data (at the resolution of the raster) for various analysis such as ordination, or incidence/frequency data for particular species can be extracted.
}
\value{Returns a list of length 2:
  \item{$grid.matrix }{A binary data.frame of species occurrences against grid cell numbers that match those in $pa.raster.}
  \item{$pa.raster }{A \code{RasterLayer} object for which species presence/absence is scored in $grid.matrix.}
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution 6, 845-852. DOI: 10.1111/2041-210X.12361
}
\author{Greg R. Guerin & Lasse Ruokolainen
}

\seealso{\code{\link{letsR::lets.presab}}
}
\examples{
#####Preparation for this example:
library(vegan)
data(mite)
data(mite.xy)
####Usage of the function:
map.pa.matrix(mite, records="site", site.coords=mite.xy)
}
