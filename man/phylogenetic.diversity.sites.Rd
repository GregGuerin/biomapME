\name{phylogenetic.diversity.sites}
\alias{phylogenetic.diversity.sites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Phylogenetic diversity for sites
%%  ~~function to do ... ~~
}
\description{Calculates (unweighted and uncorrected) Faith's phylogenetic diversity for a set of  community sample sites.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
phylogenetic.diversity.sites(species_records, phylo.tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species_records}{A data.frame with rows as sites and columns as species (presence absence matrix).
%%     ~~Describe \code{species_records} here~~
}
  \item{phylo.tree}{A phylogenetic tree of class 'phylo' containing species in the occurrence data. Must have branch lengths.
%%     ~~Describe \code{phylo.tree} here~~
}
}
\details{Unweighted and uncorrected Faith's phylogenetic diversity for community samples organised into sites. Function takes a community matrix (could be species inventories from field plots or map grid cells) and an associated tree with branch lengths. Although coded differently, this is more or less equivalent to the 'pd' function in the 'picante' package. I have pulled it out as a stand-alone function for a community matrix from the unweighted (pd) case of the biomapME::phylogenetic.endemism function, which has point record inputs and gridded map outputs.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a vector of PD scores (list with $PD as the vector) for sites in the community matrix.
}
\references{Guerin, G.R. and Lowe, A.J. (2015) Mapping phylogenetic endemism in R using georeferenced branch extents. SoftwareX 3-4, 22-26.
\cr
\cr
Faith, D. P. (1992) Conservation evaluation and phylogenetic diversity. Biological conservation, 61, 1-10.
}
\author{Greg .R. Guerin
}

\seealso{\code{\link{phylogenetic.endemism}}
\cr
\cr
\code{\link[picante]{pd}}
}
\examples{
library(vegan)
data(mite)
#generate random phylogenetic tree
mite.tree <- ape::rtree(n=ncol(mite), tip.label=colnames(mite))
mite.PD <- phylogenetic.diversity.sites(mite, phylo.tree=mite.tree)
mite.PD
}
