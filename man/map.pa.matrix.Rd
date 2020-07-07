\name{map.pa.matrix}
\alias{map.pa.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Map Presence-absence Matrix
%%  ~~function to do ... ~~
}
\description{Given georeferenced incidence data for species, generates a binary presence/absence matrix associated with grid cells of a raster.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
map.pa.matrix(species_records, records = "single", site.coords, species = "SPECIES", longitude = "LONGITUDE", latitude = "LATITUDE", frame.raster, deg.resolution = c(0.25, 0.25), extent.vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species_records}{A data.frame, either with:
  a) rows as individual species records, and columns that include fields for species name, longitude and latitude (see ‘species’, ‘longitude’, ‘latitude’ below); or
  b) rows as sites and columns as species, in which case ‘site.coords’ (below) must also be supplied.
%%     ~~Describe \code{species_records} here~~
}
  \item{records}{Are the species_records in single/long format (the default, records="single") or in site-based/short format (records="site")?
%%     ~~Describe \code{records} here~~
}
  \item{site.coords}{For site-based data (records="site"), a data.frame with rows as sites (/field plots) (names match the row names of species_records) and their geographic coordinates: first column must be x/longitude, second column y/latitude.
%%     ~~Describe \code{site.coords} here~~
}
  \item{species}{For records="single" (i.e. individual occurrence data); what colname in the supplied species_records contains species names?
%%     ~~Describe \code{species} here~~
}
  \item{longitude}{For records="single"; what colname in the supplied species_records contains longitude values?
%%     ~~Describe \code{longitude} here~~
}
  \item{latitude}{For records="single"; what colname in the supplied species_records contains latitude values?
%%     ~~Describe \code{latitude} here~~
}
  \item{frame.raster}{An existing rasterLayer object. User can elect to supply a raster, in which case presences and absences are scored for grid cells in the raster. If not specified, the function generates a raster based on default or specified extent and resolution.
%%     ~~Describe \code{frame.raster} here~~
}
  \item{deg.resolution}{Arguments specifying the map resolution (in degrees) the user wishes the calculations and mapping to use. If none are specified, default resolution (0.25) is used. If a frame.raster is specified, this argument is ignored.
%%     ~~Describe \code{deg.resolution} here~~
}
  \item{extent.vector}{Argument specifying the map extent (in degrees) the user wishes the calculations and mapping to use. If none are specified, default extent (data extent) is used. If a frame.raster is specified, this argument is ignored.
%%     ~~Describe \code{extent.vector} here~~
}
}
\details{This function generates a binary species presence/absence matrix associated with a raster layer based on georeferenced incidence data. This is a data processing step for mapping various biodiversity metrics onto raster layers. The outputs can be used as inputs into these functions, or if desired they can be used like site-based data (at the resolution of the raster) for various analysis such as ordination, or incidence/frequency data for particular species can be extracted.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a list of length 2:
  \item{$grid.matrix }{A binary data.frame of species occurrences against grid cell numbers that match those in $pa.raster.}
  \item{$pa.raster }{A ‘raster’ object for which species presence/absence is scored in $grid.matrix.}
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

\seealso{\code{\link{letsR::lets.presab}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#####Preparation for this example:
library(vegan)
data(mite)
data(mite.xy)
####Usage of the function:
map.pa.matrix(mite, records="site", site.coords=mite.xy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
