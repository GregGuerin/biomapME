\name{pe.null.test}
\alias{pe.null.test}
\title{Phylogenetic endemism null test
}
\description{Taking the outputs from the \code{phylogenetic.endemism} function, tests whether observed phylogenetic diversity/endemism is higher than expected, using non-parametric methods.
}
\usage{
pe.null.test(phylogenetic.endemism.output, nrep = 100, outlier.range = 1.5,
pe.type = "weighted", phylo.tree)
}
\arguments{
  \item{phylogenetic.endemism.output}{The returned object from a \code{\link{phylogenetic.endemism}} run, a list of length 3 (for PD, \code{pe.type="unweighted"}) or 7 (for PE: \code{pe.type="weighted"}): \code{$PD} (numeric vector of phylogenetic diversity/endemism scores for grid cells); \code{$PD_raster} (raster layer of PD/PE scores); \code{$grid.matrix} (species against cells occurrence matrix - row names match $PD names); \code{$ranges} (numeric vector of range weights for phylogenetic branches); \code{$phyloMatrix} (branches against species matrix); \code{$phylo.cell.matrix} (branches against cells matrix); \code{$edge.lengths} (numeric vector of edge.length in phylo.tree).
}
  \item{nrep}{Desired number of replicates when generating a null distribution from a random draw of species. Default is 100 for speed (slow for large datasets), but at least 1000 is recommended to ensure smooth null distributions and useful p-values.
}
  \item{outlier.range}{The multiplier of the interquartile range, outside of which scores are considered to be categorical outliers. Default is the typical value of 1.5.
}
  \item{pe.type}{Refers back to how phylogenetic.endemism was run; alternatively set to pe.type="unweighted" for phylogenetic diversity with edge length uweighted by range size if this was the original setting (in which case 'phylo.tree' must be supplied).
}
  \item{phylo.tree}{The phylogenetic tree that was used to run \code{phylogenetic.endemism}, only required for testing unweighted pd (\code{pe.type="unweighted"}).
}
}
\details{With the outputs from the \code{phylogenetic.endemism} function, performs the following tests: \strong{1)} non-parametric significance test as to whether observed phylogenetic diversity/endemism is higher or lower than expected, given species richness (and observed species frequencies); \strong{2)} Identifies and maps outliers (i.e. in terms of map grid cells that have higher or lower PD/PE) based on quantiles. As categorical: whether score lies more than 1.5 (or other user-defined amount) times outside the interquartile range; as continuous: the factor of the interquartile by which observed values differ from the median / 50 percent quantile). Returns vectors of values plus raster maps.

Raw phylogenetic diversity/endemism scores are biased both by the completeness of species sampling and species richness itself. This function does not correct the scores per se, but compares them to a null distribution. This is achieved by making replicate random draws from the species pool based on the observed species richness (i.e. same number of species) and the actual species frequencies (more frequent species more likely to be drawn). The distribution of the resulting set of null scores is compared to the observed scores and subsequently grid cells can be mapped as higher or lower than expected (based on significance testing and comparison to null quantiles).
}
\value{Returns a list with following elements. Plots of observed and expected phylogenetic diversity/endemism against species richness are generated, as are plots of the generated rasters.
  \item{$Quantile.25 }{Vector of expected lower interquartile range for a given species richness.}
  \item{$Quantile.75 }{Vector of expected upper interquartile range for a given species richness.}
  \item{$out.above.below }{Numeric vector assigning outliers categorically, with -1 for lower outlier, 0 non-outlier, 1 for upper outlier, for each grid cell.}
  \item{$out.above.below.raster }{Map of $out.above.below}
  \item{#$out.continuous }{Scores for each grid cell = by what factor of the interquartile range observed phylogenetic diversity/endemism differs from the median / 50 percent quantile.}
  \item{$out.continuous.raster }{Map of $out.continuous}
  \item{$P.above }{Vector of p-values for grid cells having higher than expected phylogenetic diversity/endemism. i.e. very low score = likely to be higher than expected; very high scores = likely to be lower than expected.}
  \item{$P.above.raster }{Map of $P.above with default 2-tailed colour scheme}
  \item{$richness }{Vector of simple observed species richness scores for cells.}
}
\references{Guerin, G.R. and Lowe, A.J. (submitted) Mapping phylogenetic endemism in R using georeferenced branch extents.
}
\author{Greg R. Guerin
}

\seealso{\code{\link{phylogenetic.endemism}}
}
\examples{
library(vegan)
data(mite)
data(mite.xy)
#for example, generate phylogenetic tree of species with random relationships and branch lengths
mite.tree <- ape::rtree(n=ncol(mite), tip.label=colnames(mite)) 
mite.PE <- phylogenetic.endemism(mite, records="site", site.coords=mite.xy, 
sep.comm.spp="none", phylo.tree=mite.tree, sep.phylo.spp="none", weight.type="cell")
pe.mite.test <- pe.null.test(mite.PE)
pe.mite.test
}
