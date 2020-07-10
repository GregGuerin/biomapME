\name{endemism.null.test}
\alias{endemism.null.test}
\title{Endemism Null Test
}
\description{Taking the outputs from the 'weighted.endemism' function, tests whether observed endemism is higher than expected, using non-parametric methods.
}
\usage{
endemism.null.test(weighted.endemism.output, nrep = 100, outlier.range = 1.5)
}
\arguments{
  \item{weighted.endemism.output}{The returned object from a 'weighted.endemism' run (or equivalent constructed separately). Consists of a list of length 4: $WE (named vector of endemism scores for grid cells), $WE_raster (raster layer of endemism scores), $weights (named vector of range weights for species used in the calculation of weighted endemism); $grid.matrix (species against cells occurrence matrix - row names must match $WE names, column names must match $weights names). If 'corrected weighted endemsim' was specified in the weighted.endemism run, this function will not accept the output. This is because correcting via dividing endemism by species richness is an alternative to this function, which uses more sophisticated methods for testing whether endemism is different than that expected for a given species richness.
}
  \item{nrep}{The desired number of replicates when generating a null distribution from a random draw of species. Default is 100 for speed, but at least 1000 is recommended to ensure smooth null distributions and useful p-values.
}
  \item{outlier.range}{Times outside the (null) interquartile range beyond which a record is considered an outlier.
}
}
\details{With the outputs from the 'weighted.endemism' function, performs the following tests:

1) Non-parametric significance test as to whether observed endemism is higher or lower than expected, given species richness (and observed species frequencies).

2) Identifies and maps outliers (i.e. in terms of map grid cells that have higher or lower endemism) based on quantiles. As categorical: whether endemism score lies more than 1.5 (or other user-defined amount) times outside the interquartile range; as continuous: the factor of the interquartile by which observed values differ from the median / 50 percent quantile). Returns vectors of values plus raster maps.

Raw weighted endemism scores are biased both by the completeness of species sampling and species richness itself. Correcting by dividing by the observed number of species ('corrected weighted endemism' of Crips et al. 2001) is a proposed correction, but the relationship between endemism scores and species richness is not linear under a null model (random species draws), as increasingly infrequent species are drawn as richness increases, thereby increasing CWE. While correcting endemism scores in a more sophisticaed way is possible, this function does not correct the scores per se, but compares them to a null distribution. This is achieved by making replicate random draws from the species pool based on the observed species richness (i.e. same number of species) and the actual species frequencies (more frequent species more likely to be drawn). The distribution of the resulting set of null endemsim scores is compared to observed endemism and subsequently grid cells can be mapped as higher or lower than expected (based on significance testing and comparison to null quantiles).
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a list with following elements. Plots of observed and expected endemism against species richness are generated, as are plots of the generated rasters.
  \item{$Quantile.25 }{Vector of expected lower interquartile range for a given species richness.}
  \item{$Quantile.75 }{Vector of expected upper interquartile range for a given species richness.}
  \item{$out.above.below }{Vector assigning outliers categorically, with -1 for lower outlier, 0 non-outlier, 1 for upper outlier for each grid cell.}
  \item{$out.above.below.raster }{Map of $out.above.below.}
  \item{$out.continuous }{Scores for each grid cell = by which factor of the interquartile range observed endemism differs from the median / 50 percent quantile.}
  \item{$out.continuous.raster }{Map of $out.continuous.}
  \item{$P.above }{Vector of p-values for grid cells having higher than expected endemism. Very low scores = more likely to be higher than expected; very high scores = more likely to be lower than expected.}
  \item{$P.above.raster }{Map of $P.above.}
  \item{$richness }{Simple observed species richness scores for cells.}
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12361
}
\author{Greg R. Guerin & Lasse Ruokolainen
}

\seealso{
  \code{\link{weighted.endemism}}
}
\examples{
library(vegan)
data(mite)
data(mite.xy)

endemism_mydata <- weighted.endemism(mite, site.coords=mite.xy, records="site")

endemism.test.example <- endemism.null.test(endemism_mydata)
}