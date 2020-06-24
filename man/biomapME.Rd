\name{biomapME}
\alias{biomapME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Biodiversity Mapping and Macroecology
%%  ~~function to do ... ~~
}
\description{The biomapME package contains a range of R functions for manipulating, analysing and rasterizing ecological data for use in macroecological analysis and mapping of biodiversity metrics. The functions are primarily geared towards generating (and analysing) gridded biodiversity maps based on various metrics from species incidence (and environmental) data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\details{biomapME is currently under active development. Additional metrics/functions/functionality and changes to the way some functions are structured are expected over time.\cr
The main purpose of the package is to provide self-contained functions in the R environment for taking simple input data (i.e. georeferenced species records) and calculating important biodiversity metrics, which are returned as automatically generated raster maps as well as in simpler numeric formats.\cr
\cr
At present, biomapME contains functions for calculating:\cr
1) Species Richness (simple and uncorrected)\cr
2) Range-rarity Richness\cr
3) Georeferenced weighted endemism\cr
4) Phylogenetic range-rarity\cr
5) Georeferenced phylogenetic endemism\cr
6) Non-parametric tests for 2):5) as to whether observed diversity is higher than expected\cr
7) Phylogenetic Diversity for a set of sample sites\cr
8) Presence-absence matrix for species occurrences within map grid cells\cr
9) Various range metrics for batch processed multiple species (number of occupied grid cells, extent of occurrence (range span, range area), longitudinal/latitudinal range)\cr
10) Batch processed habitat fragmentation/class statistics for circular buffers around focal sites\cr
11) Converting plot based species occurrence matrices into individual species records\cr
\cr
Additional functionality planned includes:\cr
1) Beta diversity metrics (with sliding scale of spatial comparison and flexible metrics)\cr
2) Incorporate existing non-parametric sampling corrections\cr
3) Extending more of the functions to suit spatial coordinate systems other than longlat (current default for most functions)
4) Increasing the range of input data formats across all function, for example species data in spatial objects

%%  ~~ If necessary, more details than the description above ~~
}
\references{Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution 6, 845–852. DOI: 10.1111/2041-210X.12361\cr
\cr
Guerin, G.R. & Lowe, A.J. (2015) Sum of inverse range-sizes ('SIR'), a biodiversity metric with many names and interpretations. Biodiversity and Conservation. DOI 10.1007/s10531-015-0977-6
%% ~put references to the literature/web site here ~
}
\author{Greg R. Guerin
%%  ~~who you are~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
