\name{convert.site.data}
\alias{convert.site.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert species records from sites with spatial coordinates to individual records (i.e. long-format).
%%  ~~function to do ... ~~
}
\description{This function converts species occurrence data from site- to individual species record-based. This
conversion could be useful in any situation where individual species records are required, but was designed primarily
to allow this data format as input to other biodiversity mapping functions, in which case this conversion makes plot
-based data compatible with the format of herbarium or musem records, for example.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
convert.site.data(species_records, site.coords)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species_records}{A data.frame with rows as sites and columns as species. Values greater than 0 are interpreted as indicating presence (abundances are not retained).}
%%     ~~Describe \code{x} here~~
 \item{site.coords}{A data.frame with the sites (/field plots) that match the column names of species_records as rows with longlat/XY coordinates.}
}
\details{This function is partially a reverse operation to simba::mama, which converts long-format (individual) species records into matrix format.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a data.frame with three columns, "SPECIES", "LONGITUDE", "LATITUDE".
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12361
%% ~put references to the literature/web site here ~
}
\author{Greg R. Guerin & Lasse Ruokolainen
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mama}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(mite)
data(mite.xy)
mite_records <- convert.site.data(mite, mite.xy)
head(mite_records)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
