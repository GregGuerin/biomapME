\name{weighted.endemism}
\alias{weighted.endemism}
\title{Range rarity richness, georeferenced weighted endemism and species richness
}
\description{Calculates (taxonomic / species) sum of inverse range-sizes ('SIR': species richness inversely weighted by species ranges) metrics across gridded maps using single or site-based point records.
}
\usage{
weighted.endemism(species_records, records = "single", site.coords,
species = "SPECIES", longitude = "LONGITUDE", latitude = "LATITUDE",
frame.raster, deg.resolution = c(0.25, 0.25), extent.vector, type = "weighted",
plot.raster = TRUE, own.weights, weight.type = "cell", geo.type = "cell",
geo.calc = "max.dist", outlier_pct = 100, verbose = TRUE, own.grid.matrix)
}
\arguments{
  \item{species_records}{A data.frame, either with:
\strong{a)} rows as individual species records, and columns that include fields for species name, longitude and latitude (see species, longitude, latitude below); or
\strong{b)} rows as sites and columns as species, in which case site.coords (below) must also be supplied.
}
  \item{records}{Are the species_records in single/long format (the default, "single") or in site-based/short format (\code{records="site"})?
}
  \item{site.coords}{For site-based data (\code{records="site"}), a data.frame with the sites (/field plots) that match the column names of species_records and their longlat coordinates.
}
  \item{species}{For \code{records="single"}; which colname in the supplied species_records contains species names?
}
  \item{longitude}{For \code{records="single"}; which colname in the supplied species_records contains longitude values?
}
  \item{latitude}{For \code{records="single"}; which colname in the supplied species_records contains latitude values?
}
  \item{frame.raster}{An existing rasterLayer object the user can optionally elect to supply as the frame for calculations and mapping. If not specified, the function will generate a raster based on default or specified extent and resolution.
}
  \item{deg.resolution}{Argument specifying the map resolution in degrees the user wishes the calculations and mapping to use. If no frame is specified, an arbitrary resolution is supplied. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster)
}
  \item{extent.vector}{Arguments specifying the map extent in degrees the user wishes the calculations and mapping to use. If no frame is specified, an arbitrary extent is supplied, that of the input point data. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster)
}
  \item{type}{Either "weighted" (default), or "corrected" (corrected weighted endemism - the 'per-species' weighted endemism as per Crisp et al. (2001). This is provided for convenience, but is not particularly recommended - instead, use the outputs of weighted.endemism in endemism.null.test function to compare to null expectations of endemism given the species richness)
}
  \item{plot.raster}{Whether or not to plot the output raster with endemism scores. Either way, the raster object is stored in the output.
}
  \item{own.weights}{An optional user-supplied numeric vector of species weights for calculating endemism. Values must have names that are a complete and exact match for the species names in species_records as each species must have a weight. This optional argument is intended mainly so that the more time consuming calculation of geographic weights can be done once and the result stored and used for subsequent re-runs of the endemism calculations.
}
  \item{weight.type}{Default is "cell" (cell-based range weights), while "geo" will calculate geographic range weights. Weight.type "richness" sets the weights to 1, which is equivalent to calculating simple species richness (note setting type="corrected" in this case will give each cell a score of 1). Argument is ignored if \code{own.weights} is supplied.
}
  \item{geo.type}{Default is "cell". This argument is only considered when weight.type is set to "geo" (above), in which case \code{geo.type="cell"} calculates geographic ranges based on map grid cell centroids. This can optionally be set to \code{geo.type="point"}, in which case the geographic range weightings are calculated based on the point locations of each species.
}
  \item{geo.calc}{Default is "max.dist". This argument is only considered when when weight.type is set to "geo" and applies to both \code{geo.type="cell"} and \code{geo.type="point"}. If \code{geo.calc="max.dist"} the default maximum geographic distance ('span') between cell centroids/point locations is calculated for a species range. If geo.calc is set to "polygon", the function weights species ranges based on the area of a MCP (minimum convex polygon) that contains all points. Further arguments to this function can be included, such as changing the default outlier_pct=95 (removes outlying locations). Additionally, the \code{geo.type="point"} option is not recommended for the "polygon" method, as it is more likely to lead to errors where nearby point locations do not allow drawing of a spanning polygon. In this case, cell-centroid based calculations ensure that multiple records are spatially separated (in different cells) and that occurrences within a single cell are returned as the area of that cell.
}
  \item{outlier_pct}{For the calculation of range span or area via convex polygons (at least five records of the species), this argument can be used to remove outliers on a percentage basis via the mcp function in package adehabitat, i.e. 95 means 5 percent most outlying points are removed. Default is 100, i.e. no outlier exclusion.
}
  \item{verbose}{
}
  \item{own.grid.matrix}{User can supply a binary matrix of species against grid cell numbers, rather than this being generated within the function. The purpose of this argument is that the step can be time consuming for large datasets, so the user can return the matrix that is returned from the function in subsequent runs with different setting to improve speed. If this is supplied, a frame.raster must also be supplied that has cell numbers which match the row.names of the own.grid.matrix.
}
}
\details{This implementation of 'weighted endemism' allows alternative calculation of weights for species ranges as well as the option of user-supplied weights. Weights can be calculated based on the frequency of occurrence in grid cells, or alternatively by the geographical size of the species range, calculated in one (span) or two (area) dimensions.
}
\value{Returns a list of length 4:
  \item{$WE (/$CWE)}{Vector of 'weighted endemism' scores.}
  \item{$WE_Raster (/$CWE)}{Raster map with endemism scores.}
  \item{$weights}{A named numeric vector of weights used to calculate endemism (equivalent to range size in metres if \code{weight.type="geo"}, range size in cells if \code{weight.type="cell"} (default), or the user supplied weights if own.weights was supplied).}
   \item{$grid.matrix}{A binary data.frame of species against grid cell numbers used in the function which is returned so that it can be re-used in subsequent runs to save time.}
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. \emph{Methods in Ecology and Evolution} \strong{6}, 845-852.
}
\author{Greg R. Guerin & Lasse Ruokolainen
}
\seealso{
  \code{\link{endemism.null.test}}
}
\examples{
library(vegan)
data(mite)
data(mite.xy)

#default, range rarity richness
endemism_mydata <- weighted.endemism(mite, site.coords = mite.xy, records="site")
endemism_mydata$WE
endemism_mydata$weights

#georeferenced weighted endemism; re-apply previous outputs
endemism_mydata2 <- weighted.endemism(mite, site.coords = mite.xy, records="site", 
weight.type="geo", own.grid.matrix = endemism_mydata$grid.matrix, 
frame.raster=endemism_mydata$WE_raster)
endemism_mydata$WE
endemism_mydata$weights


#re-apply weights from previous run
endemism_mydata3 <- weighted.endemism(mite, site.coords = mite.xy, records="site", 
own.weights = endemism_mydata2$weights)


#unweighted species richness
endemism_mydata <- weighted.endemism(mite, site.coords = mite.xy, 
records="site", weight.type="richness")
endemism_mydata$WE
endemism_mydata$weights
}
