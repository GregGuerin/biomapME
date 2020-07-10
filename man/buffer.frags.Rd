\name{buffer.frags}
\alias{buffer.frags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Batch fragmentation/class stats for circular buffers
%%  ~~function to do ... ~~
}
\description{Batch calculates habitat fragmentation indices ('fragstats') within a circular buffer zone around sites.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
buffer.frags(XY, radius, vegetation.base.raster, plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY}{The x/longitude and y/latitude coordinates of a set of sites provided in either of the following formats. If the coordinates are in columns of a data.frame, then the subsetted data.frame should be given: e.g. \code{XY = df[,c("Longitude", "Latitude")]} if there are additional columns. If the coordinates are within a SpatialPointsDataFrame, the SPDF can be given for XY (the function will locate the coordinates).
%%     ~~Describe \code{XY} here~~
}
  \item{radius}{The radius of the buffer zone around each site within which fragstats are calculated, in units of metres if the 'vegetation.base.raster' is in longlat.
%%     ~~Describe \code{radius} here~~
}
  \item{vegetation.base.raster}{A \code{RasterLayer} object containing integer classes for the presence of vegetation/habitat/land classes of 1 or more kind.
}
\item{plot}{Logical. If \code{TRUE}, the input sites and the requested buffers are overplotted on the raster.
}
}
\details{This function is a wrapper for the \code{ClassStat} function from package \code{SDMTools}, which calculates fragstats for a given raster representing a vegetation/habitat matrix. The \code{ClassStat} function is applied to a user defined circular buffer area surrounding a site/coordinate, and these are batch processed for a set of multiple sites.
}
\value{Returns a list of data frames, one for each focal location, with fragstat indices as columns and vegetation/habitat classes (read from the input raster) as rows. See example below for collapsing the list to a single data.frame object. By default, NA values within the buffer (generally representing non-habitat or non-vegetation) are returned with the value 100 for 'class'. 'Class' will otherwise be the same as the categories in the habitat raster (See \code{?ClassStat})
}
\references{Wrapper for SDMTools version of fragstats, which is based on statistics calculated by fragstats, see http://www.umass.edu/landeco/research/fragstats/fragstats.html.

Jeremy VanDerWal, Lorena Falconi, Stephanie Januchowski, Luke Shoo and Collin Storlie (2014). SDMTools: Species Distribution Modelling Tools: Tools for processing data associated with species distribution modelling exercises. R package version 1.1-221. http://CRAN.R-project.org/package=SDMTools
}
\author{Greg R. Guerin
}

\seealso{\code{\link{ClassStat}}
}
\examples{
eg.rast <- raster()
extent(eg.rast) <- c(0,10, 0, 10)
res(eg.rast) <- c(0.5,0.5)
eg.rast[] <- NA
eg.rast[sample(c(1:ncell(eg.rast)), 100)] <- 1
eg.coords <- data.frame(Longitude=c(3, 5, 7), Latitude=c(4, 6, 8))
row.names(eg.coords) <- c("SiteA", "SiteB", "SiteC")
eg.coords
frags <- suppressWarnings(buffer.frags(eg.coords, radius=150000, vegetation.base.raster=eg.rast))

#TO CONDENSE INTO SINGLE DATA.FRAME
frags <- do.call(rbind, frags) 
frags
#NOTE:class = 100 is default output for NA input raster cells.
}
