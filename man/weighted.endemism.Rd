\name{weighted.endemism}
\alias{weighted.endemism}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Range-rarity Richness, Georeferenced Weighted Endemism and Species Richness
%%  ~~function to do ... ~~
}
\description{Calculates (taxonomic / species) sum of inverse range-sizes ('SIR': species richness inversely weighted by species ranges) metrics across gridded maps using single or site-based point records.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
weighted.endemism(species_records, records = "single", site.coords, species = "SPECIES", longitude = "LONGITUDE", latitude = "LATITUDE", frame.raster, deg.resolution = c(0.25, 0.25), extent.vector, type = "weighted", plot.raster = TRUE, own.weights, weight.type = "cell", geo.type = "cell", geo.calc = "max.dist", outlier_pct = 100, verbose = TRUE, own.grid.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species_records}{A data.frame, either with:
#a) rows as individual species records, and columns that include fields for species name, longitude and latitude (see species, longitude, latitude below); or
#b) rows as sites and columns as species, in which case site.coords (below) must also be supplied.
%%     ~~Describe \code{species_records} here~~
}
  \item{records}{Are the species_records in single/long format (the default, "single") or in site-based/short format (records="site")?
%%     ~~Describe \code{records} here~~
}
  \item{site.coords}{For site-based data (records="site"), a data.frame with the sites (/field plots) that match the column names of species_records and their longlat coordinates.
%%     ~~Describe \code{site.coords} here~~
}
  \item{species}{For records="single"; what colname in the supplied species_records contains species names?
%%     ~~Describe \code{species} here~~
}
  \item{longitude}{For records="single"; what colname in the supplied species_records contains longitude values?
%%     ~~Describe \code{longitude} here~~
}
  \item{latitude}{For records="single"; what colname in the supplied species_records contains latitude values?
%%     ~~Describe \code{latitude} here~~
}
  \item{frame.raster}{An existing rasterLayer object the user can optionally elect to supply as the frame for calculations and mapping. If not specified, the function will generate a raster based on default or specified extent and resolution.
%%     ~~Describe \code{frame.raster} here~~
}
  \item{deg.resolution}{Argument specifying the map resolution in degrees the user wishes the calculations and mapping to use. If no frame is specified, an arbitrary resolution is supplied. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster)
%%     ~~Describe \code{deg.resolution} here~~
}
  \item{extent.vector}{Arguments specifying the map extent in degrees the user wishes the calculations and mapping to use. If no frame is specified, an arbitrary extent is supplied, that of the input point data. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster)
%%     ~~Describe \code{extent.vector} here~~
}
  \item{type}{Either "weighted" (default), or "corrected" (corrected weighted endemism - the 'per-species' weighted endemism as per Crisp et al. (2001). This is provided for convenience, but is not particularly recommended - insterad, use the outputs of weighted.endemism in endemism.null.test function to compare to null expectations of endemism given the species richness)
%%     ~~Describe \code{type} here~~
}
  \item{plot.raster}{Whether or not to plot the output raster with endemism scores. Either way, the raster object is stored in the output.
%%     ~~Describe \code{plot.raster} here~~
}
  \item{own.weights}{An optional user-supplied numeric vector of species weights for calculating endemism. Values must have names that are a complete and exact match for the species names in species_records as each species must have a weight. This optional argument is intended mainly so that the more time consuming calculation of geographic weights can be done once and the result stored and used for subsequent re-runs of the endemism calculations.
%%     ~~Describe \code{own.weights} here~~
}
  \item{weight.type}{Default is "cell" (cell-based range weights), while "geo" will calculate geographic range weights. Weight.type "richness" sets the weights to 1, which is equivalent to calculating simple Species Richness (note setting type="corrected" in this case will give each cell a score of 1). Argument is ignored if own.weights is supplied.
%%     ~~Describe \code{weight.type} here~~
}
  \item{geo.type}{Default is "cell". This argument is only considered when weight.type is set to "geo" (above), in which case geo.type="cell" calculates geographic ranges based on map grid cell centroids. This can optionally be set to geo.type="point", in which case the geographic range weightings are calculated based on the point locations of each species.
%%     ~~Describe \code{geo.type} here~~
}
  \item{geo.calc}{Default is "max.dist". This argument is only considered when when weight.type is set to "geo" and applies to both geo.type="cell" and geo.type="point". If geo.calc="max.dist" the default maximum geographic distance (‘span’) between cell centroids/point locations is calculated for a species range. If geo.calc is set to "polygon", the function weights species ranges based on the area of a MCP (minimum convex polygon) that contains all points. Further arguments to this function can be included, such as changing the default outlier_pct=95 (removes outlying locations). Additionally, the geo.type="point" option is not recommended for the "polygon" method, as it is more likely to lead to errors where nearby point locations do not allow drawing of a spanning polygon. In this case, cell-centroid based calculations ensure that multiple records are spatially separated (in different cells) and that occurrences within a single cell are returned as the area of that cell.
%%     ~~Describe \code{geo.calc} here~~
}
  \item{outlier_pct}{For the calculation of range span or area via convex polygons (at least 5 records of the species), this argument can be used to remove outliers on a percentage basis via the mcp function in package adehabitat, i.e. 95 means 5 percent most outlying points are removed. Default is 100, i.e. no outlier exclusion.
%%     ~~Describe \code{outlier_pct} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{own.grid.matrix}{User can supply a binary matrix of species against grid cell numbers, rather than this being generated within the function. The purpose of this argument is that the step can be time consuming for large datasets, so the user can return the matrix that is returned from the function in subsequent runs with different setting to improve speed. If this is supplied, a frame.raster must also be supplied that has cell numbers which match the row.names of the own.grid.matrix.
%%     ~~Describe \code{own.grid.matrix} here~~
}
}
\details{This implementation of 'weighted endemism' allows alternative calculation of weights for species ranges as well as the option of user-supplied weights. Weights can be calculated based on the frequency of occurrence in grid cells, or alternatively by the geographical size of the species range, calculated in one (span) or two (area) dimensions.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a list of length 4:
  \item{$WE (/$CWE)}{Vector of 'weighted endemism' scores.}
  \item{$WE_Raster (/$CWE)}{Raster map with endemism scores.}
  \item{$weights}{A named numeric vector of weights used to calculate endemism (equivalent to range size in metres if weight.type="geo", range size in cells if weight.type="cell" (default), or the user supplied weights if own.weights was supplied).}
   \item{$grid.matrix}{A binary data.frame of species against grid cell numbers used in the function which is returned so that it can be re-used in subsequent runs to save time.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution 6, 845–852. DOI: 10.1111/2041-210X.12361
%% ~put references to the literature/web site here ~
}
\author{Greg R. Guerin & Lasse Ruokolainen
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{endemism.null.test}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(vegan)
data(mite)
data(mite.xy)
endemism_mydata <- weighted.endemism(mite, site.coords = mite.xy, records="site")
endemism_mydata2 <- weighted.endemism(mite, site.coords = mite.xy, records="site", weight.type="geo", own.grid.matrix = endemism_mydata$grid.matrix, frame.raster=endemism_mydata$WE_raster)
endemism_mydata3 <- weighted.endemism(mite, site.coords = mite.xy, records="site", own.weights = endemism_mydata2$weights)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ endemism }
\keyword{ range-rarity }
\keyword{ species richness }% __ONLY ONE__ keyword per line
