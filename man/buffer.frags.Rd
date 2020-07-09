\name{buffer.frags}
\alias{buffer.frags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Batch fragmentation/class stats for circular buffers
%%  ~~function to do ... ~~
}
\description{Batch calculates habitat fragmentation indices (‘fragstats’) within a circular buffer zone around sites.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
buffer.frags(XY, radius, vegetation.base.raster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY}{The x/longitude and y/latitude coordinates of a set of sites provided in either of the following formats. If the coordinates are in columns of a data.frame, then the subsetted data.frame should be given: e.g. XY=df[,c("Longitude", "Latitude")] if there are additional columns. If the coordinates are within a SpatialPointsDataFrame, the SPDF can be given for XY (the function will locate the coordinates).
%%     ~~Describe \code{XY} here~~
}
  \item{radius}{The radius of the buffer zone around each site within which fragstats are calculated, in units of metres if the 'vegetation.base.raster' is in longlat.
%%     ~~Describe \code{radius} here~~
}
  \item{vegetation.base.raster}{A \code{RasterLayer} object containing integer classes for the presence of vegetation/habitat/land classes of 1 or more kind.
%%     ~~Describe \code{vegetation.base.raster} here~~
}
}
\details{This function is a wrapper for the \code{ClassStat} function from package \code{SDMTools}, which calculates fragstats for a given raster representing a vegetation/habitat matrix. The ClassStat function is applied to a user defined circular area surrounding a site/coordinate, and these are batch processed for multiple sites. The circular buffer is generated through the function \code{raster::buffer}, a process that can be slow for high resolution rasters or large areas.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a list of data frames, one for each focal location, with fragstat indices as columns and vegetation/habitat classes (read from the input raster) as rows. See example below for collapsing the list to a single data.frame object. By default, NA values within the buffer (generally representing non-habitat or non-vegetation) are returned with the value 100 for 'class'. 'Class' will otherwise be the same as the categories in the habitat raster (See \code{?ClassStat})
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Wrapper for SDMTools version of fragstats, which is based on statistics calculated by fragstats, see http://www.umass.edu/landeco/research/fragstats/fragstats.html.

Jeremy VanDerWal, Lorena Falconi, Stephanie Januchowski, Luke Shoo and Collin Storlie (2014). SDMTools: Species Distribution Modelling Tools: Tools for processing data associated with species distribution modelling exercises. R package version 1.1-221. http://CRAN.R-project.org/package=SDMTools
%% ~put references to the literature/web site here ~
}
\author{Greg R. Guerin
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{ClassStat}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
eg.rast <- raster()
extent(eg.rast) <- c(0,10, 0, 10)
res(eg.rast) <- c(0.5,0.5)
eg.rast[] <- NA
eg.rast[sample(c(1:ncell(eg.rast)), 100)] <- 1
plot(eg.rast)
eg.coords <- data.frame(Longitude=c(3, 7), Latitude=c(4, 8))
row.names(eg.coords) <- c("SiteA", "SiteB")
eg.coords
points(eg.coords, pch=20, cex=2)
frags <- buffer.frags(eg.coords, radius=150000, vegetation.base.raster=eg.rast)

frags <- do.call(rbind, frags) #TO CONDENSE INTO SINGLE DATA.FRAME
frags #NOTE:class =100 in 'frags' is the default for NAs in the input raster (i.e. no habitat/vegetation).
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
