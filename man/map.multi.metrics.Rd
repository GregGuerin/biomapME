\name{map.multi.metrics}
\alias{map.multi.metrics}
\title{Alpha and multi-site beta species diversity within grid cells
}
\description{Given plot-based incidence data for species (i.e., presence/absence), calculates either alpha diversity (observed species richness plus a set of estimates using nonparametric methods), multi-site beta diversity associated with grid cells of a raster, based on multiple plot samples per cell, or phylogenetic diversity (observed plues estimated based on completeness of species sampling).
}
\usage{
map.multi.metrics(species_records, site.coords, alpha=TRUE, beta=TRUE, 
phylogenetic=FALSE, phylo.tree, frame.raster, deg.resolution=c(0.25,0.25), 
extent.vector, plot.raster=TRUE)
}
\arguments{
  \item{species_records}{A data.frame of species occurrences in plots with rows as sites and columns as species. Cannot make use of presence-only records.
}
  \item{site.coords}{A data.frame with rows as sites (/field plots) (names match the row names of species_records) and their geographic coordinates: first column must be x/longitude, second column y/latitude.
}
  \item{alpha}{Logical with default TRUE. Whether to calculate alpha diversity  for grid cells using nonparametric estimators.
}
  \item{beta}{Logical with default TRUE. Whether to calculate multi-site beta diversity for grid cells.
}
  \item{phylogenetic}{Logical with default FALSE. Whether to calculate phylogenetic diversity for grid cells with estimates based on completeness of species sampling.
}
  \item{phylo.tree}{A phylogenetic tree of class 'phylo' containing species in the occurrence data. Must have branch lengths.
}
  \item{frame.raster}{Optional. An existing \code{rasterLayer} object. User can elect to supply a raster. If not specified, the function generates a raster based on default or specified extent and resolution. Alpha and beta metrics are calculated based on plots in \code{species_records} that fall within grid cells of the \code{frame.raster}. 
}
  \item{deg.resolution}{Arguments specifying the map resolution (in degrees) the user wishes the calculations and mapping to use. If none are specified, default resolution (0.25) is used. If a frame.raster is specified, this argument is ignored.
}
  \item{extent.vector}{Argument specifying the map extent (in degrees) the user wishes the calculations and mapping to use. If none are specified, default extent (data extent) is used. If a frame.raster is specified, this argument is ignored.
}
  \item{plot.raster}{Logical with default TRUE. Whether to plot out generated rasters of alpha and beta diversity. Either way, results tables and rasters are returned.}
  
}
\details{This function calculates diversity among plots located within the same grid cells, for metrics that are based on multiple community observations, rather than pooled observations over the entire grid cell. \code{map.multi.metrics} uses plot-based function in packages BAT and vegan.

Alpha diversity estimation uses the \code{\link[vegan]{specpool}} function, which provides observed species richness as well as estimates of richness based on chao, jackknife and bootstrap. All estimates are returned in the table but only chao is rasterised. This complements calculation of observed unweighted species richness (see \code{\link{weighted.endemism}} for plot-based datasets.)

For multi-site beta diversity, the \code{\link[BAT]{beta.multi}} function is used to calculate diversity among the plots. Total beta diversity as well as its replacement and richness components (Cardoso et al. 2014a) are returned in the table but only total beta is rasterised.

Phylogenetic diversity estimation is based on the phylogenetic diversity of the entire species assemblage at grid cell level (i.e., all species in plots that fall in that cell) corrected by the completeness of species sampling (Cardoso et al. 2014b). In this case, PD is divided by observed/estimated species richness (based on the above alpha diversity estimation), so that less complete sampling results in a larger boost to estimated PD. This complements calculation of observed PD in grid cells (see \code{\link{phylogenetic.endemism}}) for plot-based datasets and individual plot samples (see \code{\link{phylogenetic.diversity.sites}}).

The function is suitable for mapping grid-cell level estimates of multi-site alpha, beta and phylogenetic diversity from distributed networks of vegetation plots.
}
\value{Returns a list of length 2, 4 or 6 depending on whether alpha, beta and phylogenetic are set to TRUE:
  \item{$alpha_result}{A data.frame of alpha diversity metrics against rows which are named by grid cell numbers that match those in the raster.}
  \item{$alpha_raster}{A \code{RasterLayer} object with values that represent chao estimates of species richness.}
  \item{$beta_result}{A data.frame of beta diversity metrics against rows which are named by grid cell numbers that match those in the raster. Btotal is made up of the other two components, Brepl and Brich, which represent the richness difference and species replacement componenets of beta diversity.}
  \item{$beta_raster}{A \code{RasterLayer} object with values that represent total multi-site beta diversity.}
    \item{$phylo_result}{A data.frame of phylogenetic diversity (observed and estimated) against rows which are named by grid cell numbers that match those in the raster.}
  \item{$phylo_raster}{A \code{RasterLayer} object with values that represent estimated phylogenetic diversity.}
}
\references{Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014a) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. \emph{Journal of Biogeography} \strong{41}, 749-761.

Cardoso, P., Rigal, F., Borges, P.A. & Carvalho, J.C. (2014b) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. \emph{Methods in Ecology and Evolution} \strong{5}, 452-461.

Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics} \strong{43}, 783-791.
}
\author{Greg R. Guerin
}

\seealso{\code{\link{phylogenetic.diversity.sites}
\cr
\cr
\code{\link{phylogenetic.endemism}}
\cr
\cr
\code{\link{weighted.endemism}}
\cr
\cr
\code{\link[vegan]{specpool}}
\cr\cr
\code{\link[BAT]{beta.multi}}
}
\examples{
#####Preparation for this example:
library(vegan)
data(mite)
data(mite.xy)
#generate random phylogenetic tree
mite.tree <- ape::rtree(n=ncol(mite), tip.label=colnames(mite))
####Usage of the function:
map.multi.metrics(mite, mite.xy, phylogenetic=TRUE, phylo.tree=mite.tree, deg.resolution=c(1, 1))
}
