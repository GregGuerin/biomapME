\name{biomapME}
\alias{biomapME}
\docType{package}
\title{Biodiversity Mapping and Macroecology
}
\description{The biomapME package contains a range of R functions for manipulating, analysing and rasterising ecological data for use in macroecological analysis and mapping of biodiversity metrics. The functions are primarily geared towards generating (and analysing) gridded biodiversity maps based on various metrics from species incidence (and environmental) data.
}
\details{\tabular{ll}{
Package: \tab biomapME\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2022-05-31\cr
License: \tab GPL-3\cr
}
biomapME is currently under active development. Additional metrics/functions/functionality are expected over time.\cr

The main purpose of the package is to provide self-contained functions in the R environment for taking simple input data (i.e., georeferenced presence or plot-based species records) and calculating important biodiversity metrics, which are returned as automatically generated raster maps or in numeric formats.\cr
\cr

\emph{Functions in the ausplotsR package:}

At present, biomapME contains functions for calculating:\cr

\code{\link{buffer.frags}} : Batch-processed habitat fragmentation/class statistics for circular buffers around focal sites.

\code{\link{convert.site.data}} : Converting plot based species occurrence matrices into individual species records.

\code{\link{endemism.null.test}} : Non-parametric tests as to whether observed endemism is higher than expected.

\code{\link{map.multi.metrics}} : Metrics based on sets of multiple presence/absence plots within map grid cells: 1. Nonparametric estimators of species richness. 2. Multi-site beta diversity. 3. Estimated phylogenetic diversity based on completeness of species sampling. 

\code{\link{map.pa.matrix}} : Presence-absence matrices for species occurrences within map grid cells.

\code{\link{pe.null.test}} : Non-parametric tests as to whether observed phylogenetic diversity/endemism is higher than expected.

\code{\link{phylogenetic.diversity.sites}} : Phylogenetic diversity for a set of sample sites.

\code{\link{phylogenetic.endemism}} : Phylogenetic diversity, Phylogenetic range rarity and Georeferenced phylogenetic endemism.

\code{\link{spp.range.metrics}} : Various range metrics for batch-processed multiple species (i.e, number of occupied grid cells, extent of occurrence (range span, range area), longitudinal/latitudinal range).

\code{\link{weighted.endemism}} : Species richness, Range rarity richness and Georeferenced weighted endemism.
\cr
}

\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. \emph{Methods in Ecology and Evolution} \strong{6}, 845-852.\cr
\cr
Guerin, G.R. & Lowe, A.J. (2015) Sum of inverse range-sizes ('SIR'), a biodiversity metric with many names and interpretations. \emph{Biodiversity and Conservation} \strong{24}, 2877-2882.
}
\author{Greg R. Guerin
}
