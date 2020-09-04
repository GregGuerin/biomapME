\name{spp.range.metrics}
\alias{spp.range.metrics}
\title{Species range metrics
}
\description{Calculates selected species range metrics, given a set of georeferenced locations.  Metrics include grid cell and point location based measures of Area of Occupancy and Extent of Occurrence.
}
\usage{
spp.range.metrics(species_records, species = "SPECIES", longitude = "LONGITUDE", 
latitude = "LATITUDE", coord.type = "longlat", weight.type = "cell",
geo.calc = "max.dist", outlier_pct = 100, verbose = TRUE, frame.raster,
deg.resolution = c(0.25, 0.25), extent.vector, plot.out = TRUE)
}
\arguments{
  \item{species_records}{A data.frame with rows as individual species records, and columns that include fields for species name, longitude/x and latitude/y (see species, longitude, latitude below).
}
  \item{species}{What colname in the supplied species_records contains species names?
}
  \item{longitude}{What colname in the supplied species_records contains longitude/Y values?
}
  \item{latitude}{What colname in the supplied species_records contains latitude/Y values?
}
  \item{coord.type}{Either the default of "longlat" or "custom" for any other coordinate systems. If the former, distances will be calculated in metres, if "custom", in the the same units as the XY coordinates.
}
  \item{weight.type}{Default is "cell" (cell-based range metrics), while "geo" will calculate geographic range weights (no rasters).
}
  \item{geo.calc}{Default is "max.dist". This argument is only considered when when weight.type is set to "geo". If \code{geo.calc="max.dist"} the default maximum geographic distance ('span') between point locations is calculated for a species range. If geo.calc is set to "polygon", the function weights species ranges based on the area of a MCP (minimum convex polygon) that contains all points. Further arguments to this function can be included, such as changing the default outlier_pct=95 (removes outlying locations). To calculate longitudinal/latitudinal range spans, set to: \code{geo.calc="LONG"} or \code{geo.calc="LAT"}.
}
  \item{outlier_pct}{For the calculation of range span or area via convex polygons (at least 5 records of the species), this argument can be used to remove outliers on a percentage basis via the mcp function in package adehabitat, i.e. 95 means 5 percent most outlying points are removed (the default).
}
  \item{verbose}{If TRUE (default) detail of progress by species is reported while the funciton is still running. Alternatively, set to FALSE.
}
  \item{frame.raster}{An existing rasterLayer object the user can optionally elect to supply as the frame for calculations and mapping. If not specified, the function will generate a raster based on default or specified extent and resolution (if needed).
}
  \item{deg.resolution}{Argument specifying the map resolution (in degrees for longlat) the user wishes the calculations and mapping to use. If no frame is specified, an arbitrary resolution is supplied. If a frame.raster is specified, argument is ignored (function bases mapping on the supplied raster).
}
  \item{extent.vector}{Argument specifying the map  extent (in degrees for longlat) the user wishes the calculations and mapping to use. If no frame is specified, an arbitrary extent is supplied. If a frame.raster is specified, argument is ignored (function bases mapping on the supplied raster).
}
  \item{plot.out}{TRUE or FALSE. Whether or not to plot locations and rasters for cell count range scores (\code{weight.type="cell"}). Either way, the rasters are stored in the output.
}
}
\details{
Given a list of georeferenced (longlat) species records, calculates range metrics according to user choice of number of occupied map grid cells (from a supplied raster or one automatically generated within the function), maximum span across the range, range area (area of polygon defined by occurrences) or the longitudinal/latitudinal extent (defaulting to 1 for single records).

Georeferenced distance and area range metrics can be calculated longlats (default and preferred) or any coordinate set. For non-longlat coordinates, distance and area calculations assume a flat surface (longlat calcs assume curved surface). This would be appropriate for smaller areas where distance measures won't be as biased by curvature of the Earth's surface, so that the function can be used directly on data sets where species coordinates are not recorded in longlat, e.g. fine-scale locations within plots, or UTM coordinates within one map zone etc.
}
\value{Returns a list. For cell-based frequency, contains a vector of range scores (number of cells) and a RasterStack containing a RasterLayer showing cell occupancy for each unique species. For georeferenced calculations of range, only a vector of range scores (in m or m^2 or degrees of longitude/latitude if longlat coordinates used)
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. \emph{Methods in Ecology and Evolution} \strong{6}, 845-852.
}
\author{Greg R. Guerin
}

\seealso{\code{\link{weighted.endemism}}
\cr
\cr
\code{\link[letsR]{lets.rangesize}}
}
\examples{
library(vegan)
data(mite)
data(mite.xy)
mite_records <- convert.site.data(mite, mite.xy)
spp.range.metrics(mite_records, coord.type="custom", weight.type="geo")
}
