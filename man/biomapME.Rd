\name{biomapME}
\alias{biomapME}
\title{Biodiversity Mapping and Macroecology
}
\description{The biomapME package contains a range of R functions for manipulating, analysing and rasterising ecological data for use in macroecological analysis and mapping of biodiversity metrics. The functions are primarily geared towards generating (and analysing) gridded biodiversity maps based on various metrics from species incidence (and environmental) data.
}
\details{biomapME is currently under active development. Additional metrics/functions/functionality and changes to the way some functions are structured are expected over time.\cr
The main purpose of the package is to provide self-contained functions in the R environment for taking simple input data (i.e., georeferenced presence or plot-based species records) and calculating important biodiversity metrics, which are returned as automatically generated raster maps or in numeric formats.\cr
\cr
At present, biomapME contains functions for calculating:\cr
1) Species Richness (simple and uncorrected)\cr
2) Range rarity richness\cr
3) Georeferenced weighted endemism\cr
4) Phylogenetic range rarity\cr
5) Georeferenced phylogenetic endemism\cr
6) Non-parametric tests for 2):5) as to whether observed diversity is higher than expected\cr
7) Phylogenetic diversity for a set of sample sites\cr
8) Presence-absence matrices for species occurrences within map grid cells\cr
9) Various range metrics for batch-processed multiple species (i.e, number of occupied grid cells, extent of occurrence (range span, range area), longitudinal/latitudinal range)\cr
10) Batch-processed habitat fragmentation/class statistics for circular buffers around focal sites\cr
11) Converting plot based species occurrence matrices into individual species records\cr
\cr
Additional functionality planned includes:\cr
1) Beta diversity metrics (with sliding scale of spatial comparison and flexible metrics)\cr
2) Incorporate existing non-parametric sampling corrections or tests\cr
3) Extending more of the functions to suit spatial coordinate systems other than longlat (current default for most functions)\cr
4) Increasing the range of input data formats across all function, for example species data in spatial objects
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. \emph{Methods in Ecology and Evolution} \strong{6}, 845-852.\cr
\cr
Guerin, G.R. & Lowe, A.J. (2015) Sum of inverse range-sizes ('SIR'), a biodiversity metric with many names and interpretations. \emph{Biodiversity and Conservation} \strong{24}, 2877-2882.
}
\author{Greg R. Guerin
}
