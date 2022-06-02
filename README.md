# biomapME — Biodiversity Mapping and Macroecology
R functions for generating gridded biodiversity maps from point (incidence) data


biomapME contains R functions for manipulating, analysing and rasterising ecological data for use in macroecological analysis and mapping of biodiversity metrics. The functions are primarily geared towards generating gridded biodiversity maps based on various metrics from species incidence data. A second focus is generating biodiversity metrics for plot-based community samples.

biomapME remains a package in development to add functionality. The 'original' v1.0 functions (cited in various publications) can be found in the archived version cited below. However, some aspects do not work in the archived scripts due to changes in underlying package dependencies. These glitches have been corrected in v2.0.

The main purpose of the package is to provide R functionality for taking simple input data (i.e., georeferenced species records, whether presence-only or plot-based) and calculating biodiversity metrics for raster maps (or sites). Results are generally returned as automatically generated raster maps and in numeric formats.

biomapME contains functions for calculating species richness, phylogenetic diversity, 'weighted endemism' (and its various guises including phylogenetic endemism and georeferenced range span/area implementations), with non-parametric tests for higher diversity than expected by chance. Supporting functions include calculation of various range metrics (batch-processed for multiple species and includes number of occupied grid cells, extent of occurrence (range span, range area), and longitudinal/latitudinal range), generation of species presence-absence matrices for species within map grid cells, conversion of plot based species occurrence matrices into individual species records, and batch-processed habitat fragmentation/class statistics for circular buffers around focal sites (useful for testing landscape influence on site-based metrics, for example).


The current version of biomapME imports functions from other packages, which I recommend to install first: raster, geosphere, simba, SDMTools, ape, pracma, maps, vegan, adehabitatHR, dismo, sp, BAT and PhyloMeasures (knitr,
rmarkdown, markdown suggested for vignette building). biomapME can be installed and accessed as follows:
```
library(devtools)
install_github("GregGuerin/biomapME", build_vignettes = TRUE, dependencies = FALSE)
library(biomapME)
help(biomapME)
``` 


biomapME v2.1 can be cited as:
```
Guerin, G.R. (2021) biomapME: Biodiversity Mapping and Macroecology. R package v2.1. https://github.com/GregGuerin/biomapME
```


biomapME v1.0 (archived) can be cited as:
```
Guerin, G.R. (2015) biomap: Biodiversity Mapping Package: Initial release. R package v1.0. (https://github.com/GregGuerin/biomap) (https://doi.org/10.5281/zenodo.23116)
```

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.23116.svg)](http://dx.doi.org/10.5281/zenodo.23116)

