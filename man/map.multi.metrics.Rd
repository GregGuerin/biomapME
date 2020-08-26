\name{map.multi.metrics}
\alias{map.multi.metrics}
\title{Alpha and multi-site beta species diversity within grid cells
}
\description{Given plot-based incidence data for species (i.e., presence/absence), calculates alpha diversity (observed species richness plus a set of estimates using nonparametric methods) and multi-site beta diversity associated with grid cells of a raster, based on multiple plot samples per cell.
}
\usage{
map.multi.metrics(species_records, site.coords, alpha=TRUE, beta=TRUE, 
frame.raster, deg.resolution=c(0.25,0.25), extent.vector, plot.raster=TRUE)
}
\arguments{
  \item{species_records}{A data.frame of species occurrences in plots withrows as sites and columns as species.
}
  \item{site.coords}{A data.frame with rows as sites (/field plots) (names match the row names of species_records) and their geographic coordinates: first column must be x/longitude, second column y/latitude.
}
  \item{alpha}{Logical with default TRUE. Whether to calculate alpha diversity  for grid cells using nonparametric estimators.
}
  \item{beta}{Logical with default TRUE. Whether to calculate multi-site beta diversity for grid cells.
}
  \item{frame.raster}{Optional. An existing \code{rasterLayer} object. User can elect to supply a raster. If not specified, the function generates a raster based on default or specified extent and resolution. Alpha and beta metrics are calculated based on plots in \code{species_records} that fall within grid cells of the \code{frame.raster}. 
}
  \item{deg.resolution}{Arguments specifying the map resolution (in degrees) the user wishes the calculations and mapping to use. If none are specified, default resolution (0.25) is used. If a frame.raster is specified, this argument is ignored.
}
  \item{extent.vector}{Argument specifying the map extent (in degrees) the user wishes the calculations and mapping to use. If none are specified, default extent (data extent) is used. If a frame.raster is specified, this argument is ignored.
}
  \item{plot.raster}{Logical with default TRUE. Whether to plot out generated rasters of alpha and beta diversity. Either way, results tables and rasters are returned.}
  
}
\details{This function calculates diversity among plots located within the same grid cells. It is effectively a wrapper for plot-based function in packages BAT and vegan. Alpha diversity estimation uses the \code{\link[vegan]{specpool}} function, which provides observed species richness as well as estimates of richness based on chao, jackknife and bootstrap. All estimates are returned in the table but only chao is rasterised. For multi-site beta diversity, the \code{\link[BAT]{beta.multi}} function is used to calculate diversity among the plots. Total beta diversity as well as its replacement and richness components are returned in the table but only total beta is rasterised.

The function is suitable for mapping grid-cell level estimates of alpha and beta diversity from distributed networks of vegetation plots.
}
\value{Returns a list of length 4 if both alpha and beta diversity are set to TRUE:
  \item{$alpha_result }{A data.frame of alpha diversity metrics against rows which are named by grid cell numbers that match those in the raster.}
  \item{$alpha_raster }{A \code{RasterLayer} object with values that represent chao estimates of species richness.}
  \item{$beta_result }{A data.frame of beta diversity metrics against rows which are named by grid cell numbers that match those in the raster. Btotal is made up of the other two components, Brepl and Brich, which represent the richness difference and species replacement componenets of beta diversity.}
  \item{$beta_raster }{A \code{RasterLayer} object with values that represent total multi-site beta diversity.}
}
\references{Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. \emph{Journal of Biogeography} \strong{41}, 749-761.

Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics} \strong{43}, 783-791.
}
\author{Greg R. Guerin
}

\seealso{\code{\link[vegan]{specpool}}
\cr\cr
\code{\link[BAT]{beta.multi}}
}
\examples{
#####Preparation for this example:
library(vegan)
data(mite)
data(mite.xy)
####Usage of the function:
map.multi.metrics(mite, mite.xy, deg.resolution=c(1, 1))
}
