\name{phylogenetic.endemism}
\alias{phylogenetic.endemism}
\title{Phylogenetic diversity, phylogenetic range rarity & georeferenced phylogenetic endemism
}
\description{Calculates phylogenetic endemism (Faith's phylogenetic diversity inversely weighted by the spatial range of particular branch lengths) or alternatively (unweighted) phylogenetic diversity across gridded maps using individual or site-based point records.
}
\usage{
phylogenetic.endemism(species_records, records = "single", site.coords,
species = "SPECIES", longitude = "LONGITUDE", latitude = "LATITUDE",
sep.comm.spp = " ", phylo.tree, sep.phylo.spp = "_", frame.raster,
deg.resolution = c(0.25, 0.25), extent.vector, pe.type = "weighted",
plot.raster = TRUE, weight.type = "cell", outlier_pct = 100, own.weights,
own.grid.matrix, own.phylo.cell.matrix, own.phyloMatrix, pd.standard = FALSE)
}
\arguments{
  \item{species_records}{A data.frame, either with: a) rows as individual species records, and columns that include fields for species name, longitude and latitude (see 'species', 'longitude', 'latitude' below); or b) rows as sites and columns as species, in which case 'site.coords' (below) must also be supplied.
}
  \item{records}{Are the species_records in single/long format (the default, \code{records="single"}) or in site-based/short format (\code{records="site"}).
}
  \item{site.coords}{For site-based data (records="site"), a data.frame with rows as sites (/field plots) (names match the. row names of species_records) and their geographic (longlat) coordinates: first column must be x/longitude, second column y/latitude
}
  \item{species}{For \code{records="single"} (i.e. individual occurrence data); what colname in the supplied species_records contains species names?
}
  \item{longitude}{for \code{records="single"}; which colname in the supplied species_records contains longitude values?
}
  \item{latitude}{For \code{records="single"}; which colname in the supplied species_records contains latitude values?
}
  \item{sep.comm.spp}{The genus_species separating character in the community data. If there is none (i.e. taxon names represented by a single 'word'), argument should be set to \code{sep.comm.spp="none"}. The purpose of this argument is that it is common for separators to differ between community data and phylogenetic trees due to the different processing functions applied, and to avoid having to reformat the data, this argument allows them to be matched within the function. Default is a space.
}
  \item{phylo.tree}{A phylogenetic tree of class 'phylo' containing species in the occurrence data. Must have branch lengths.
}
  \item{sep.phylo.spp}{The genus_species separating character in the phylo.tree (tip.labels). If there is none (i.e. taxon names represented by a single 'word'), argument should be set to \code{sep.phylo.spp="none"}. See sep.comm.spp argument for purpose. Default for phylo.tree is an underscore, as this is a common format in community ecology using e.g. phylomatic trees.
}
  \item{frame.raster}{An existing rasterLayer object. User can elect to supply a raster as the frame for calculations and mapping. If not specified, the function will generate a raster based on default or specified extent and resolution.
}
  \item{deg.resolution}{Argument specifying the map resolution in degrees the user wishes the calculations and mapping to use. If none are specified, default resolution (0.25) is used. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster).
}
  \item{extent.vector}{Arguments specifying the map extent in degrees the user wishes the calculations and mapping to use. If none are specified, default extent (data extent) is used. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster).
}
  \item{pe.type}{Either "weighted" (default; = phylogenetic endemism), or "unweighted" (= phylogenetic diversity).
}
  \item{plot.raster}{Logical. Whether or not to plot the output raster of endemism scores. Whether plotted or not, a raster object is stored in the returned output.
}
  \item{weight.type}{Range wieght type for calculation of phylogenetic endemism. Default is "cell" (grid cell-based range weights:: phylogenetic range rarity), while "geo" will calculate georeferenced range 'span' weights (georeferenced phylogenetic endemism). Argument is ignored if \code{pe.type="unweighted"} or if \code{own.weights} are provided.
}
  \item{outlier_pct}{For the calculation of range 'span' via convex polygons, this argument can be used to remove outliers on a percentage basis via the 'mcp' function in package adehabitatHR, i.e. 95 means 5 percent most outlying points are removed. By default, 100, i.e., outliers are not excluded.
}
  \item{own.weights}{An optional user-supplied numeric vector of branch length weights for calculating phylogenetic endemism. This argument is intended mainly so that time consuming calculation of georeferenced weights can be done once and the result stored and used for subsequent re-runs. Alternatively, user-calculated weights can be used in subsequent runs. The weights must be matched to a matrix representation of phylo.tree (see below).
}
  \item{own.grid.matrix}{Optional user supplied binary matrix of species against grid cell numbers, rather than this being generated within the function. The purpose of this argument is to save time for repeat runs, given the step can be time consuming for large datasets. Must be provided if \code{own.phylo.cell.matrix} is supplied-- although the function no longer requires the grid.matrix, it is still needed for inclusion in the outputs, so that they can be subsequently fed into the pe.null.test function (which does requires the grid.matrix...). If supplied, a frame.raster must also be supplied that has cell numbers which match the row.names of the own.grid.matrix. Assumes that genus/species words are separated by the sep.phylo.spp, as this is how the grid.matrix is returned from the function originally.
}
  \item{own.phylo.cell.matrix}{User can supply the matrix representing occurrences of phylogenetic branches in map grid cells from previous runs to save computation time (i.e. \code{$phylo.cell.matrix} of a previous \code{phylogenetic.endemism output}). A matching 'frame.raster' must also be supplied.
}
  \item{own.phyloMatrix}{If user provides the above \code{own.phylo.cell.matrix} in a repeat run of the function, this must also be supplied, as the function still defaults to returning 'phylo.matrix' even though the function now skips its generation. e.g. \code{$phyloMatrix} of an previous phylogenetic.endemism.R run).
}
  \item{pd.standard}{Logical with default FALSE. Applies only to the calculation of unweighted phylogenetic diversity (i.e., \code{pe.type="unweighted"}). If set to TRUE, PD is standardised by species richness. See Details.
}
}
\details{This implementation of phylogenetic endemism for incidence data within map grid cells allows alternative calculation of weights for branch length geographic ranges. Weights can be calculated based on the frequency of occurrence in grid cells, or alternatively by the georeferenced span of the range.

Unweighted phylogenetic diversity can also be selected. In this case, if \code{pd.standard = TRUE}, the function becomes a rasterised wrapper for the \code{\link[PhyloMeasures]{pd.query}} function in package PhyloMeasures, which uses a null model to standardise PD for species richness (given PD is additive by species) by subtracting the mean and dividing by the standard deviation of expected PD. 
}

\value{Returns a list of length 3 (for \code{pe.type="unweighted"} i.e. phylogenetic diversity/PD) or length 7 (for \code{pe.type="weighted"} i.e. phylogenetic endemism/PE). For both PD and PE:
  \item{$PD }{Vector of phylogenetic diversity/endemism scores. }
  \item{$PD_raster }{Raster map with phylogenetic diversity/endemism scores.}
    \item{$grid.matrix }{A binary data.frame of species against grid cell numbers used in the function which is returned so that it can be re-used to save computation time, and because it is required in the downstream pe.null.test.R function to calculated species sample probabilities. Note that this only includes species that are found in the phylo.tree.}
For PE only:
  \item{$ranges }{A named numeric vector of weights used to calculate endemism (equivalent to range size in metres if \code{weight.type="geo"}, range size in cells if weight.type="cell" (default), or the user supplied weights if own.weights was supplied) (skipped if \code{pe.type="unweighted"}). }
  \item{$phyloMatrix }{A binary matrix representation of the phylo.tree, returned mainly for convenience for downstream use (skipped if \code{pe.type="unweighted"}). }
  \item{$phylo.cell.matrix }{A binary matrix recording the presence of particular phylogenetic branches in map grid cells used in the function, which is returned for re-use in subsequent runs for efficiency (e.g. with different weights), (skipped if pe.type="unweighted"). }
  \item{$edge.lengths }{A numeric vector of edge lengths from the phylo.tree }
}
\references{Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. \emph{Biological Conservation} \strong{61}, 1-10.
\cr
\cr
Guerin, G.R. & Lowe, A.J. (2015) Mapping phylogenetic endemism in R using georeferenced branch extents. \emph{SoftwareX} \strong{3-4}, 22-26.
\cr
\cr
Rosauer, D., Laffan, S.W., Crisp, M.D., Donnellan, S.C. & Cook, L.G. (2009) Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. \emph{Molecular Ecology} \strong{18}, 4061-4072.
\cr
\cr
Tsirogiannis, C. & B. Sandel. (2016) PhyloMeasures: A package for computing phylogenetic biodiversity measures and their statistical moments. \emph{Ecography} \strong{39}, 709-714,
}
\author{Greg R. Guerin
}

\seealso{\code{\link{pe.null.test}}
\cr
\cr
\code{\link{phylogenetic.diversity.sites}}
\cr
\cr
\code{\link{map.multi.metrics}}
}
\examples{
library(vegan)
data(mite)
data(mite.xy)

#generate random phylogenetic tree of species in dataset
mite.tree <- ape::rtree(n=ncol(mite), tip.label=colnames(mite)) 

#calculate phylogenetic range rarity
mite.PE <- phylogenetic.endemism(mite, records="site", site.coords=mite.xy, 
sep.comm.spp="none", phylo.tree=mite.tree, sep.phylo.spp="none", weight.type="cell")
mite.PE$PD

#calculate richness-standardised phylogenetic diversity
mite.PD <- phylogenetic.endemism(mite, records="site", site.coords=mite.xy, 
sep.comm.spp="none", phylo.tree=mite.tree, sep.phylo.spp="none", 
pe.type="unweighted", pd.standard=TRUE)
mite.PD$PD
}
