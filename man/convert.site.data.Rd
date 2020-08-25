\name{convert.site.data}
\alias{convert.site.data}
\title{Convert species records from sites with spatial coordinates to individual records (i.e. long-format).
}
\description{This function converts species occurrence data from site- to individual species record-based. This
conversion could be useful in any situation where individual species records are required, but was designed primarily
to allow this data format as input to other biodiversity mapping functions, in which case this conversion makes plot
-based data compatible with the format of herbarium or musem records, for example.
}
\usage{
convert.site.data(species_records, site.coords)
}
\arguments{
  \item{species_records}{A data.frame with rows as sites and columns as species. Values greater than 0 are interpreted as indicating presence (abundances are not retained).}
  \item{site.coords}{A data.frame with the sites (/field plots) that match the column names of species_records as rows with longlat/XY coordinates.}
}
\details{This function is partially a reverse operation to simba::mama, which converts long-format (individual) species records into matrix format.
}
\value{Returns a data.frame with three columns, "SPECIES", "LONGITUDE", "LATITUDE".
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. \emph{Methods in Ecology and Evolution} \strong{6}, 845-852.
}
\author{Greg R. Guerin & Lasse Ruokolainen
}
\seealso{\code{\link{mama}}
}
\examples{
library(vegan)
data(mite)
data(mite.xy)
mite_records <- convert.site.data(mite, mite.xy)
head(mite_records)
}