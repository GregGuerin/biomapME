pe.null.test <- function(phylogenetic.endemism.output, nrep = 100, outlier.range = 1.5, pe.type="weighted", phylo.tree)
#
#Description --
#
#Taking the outputs from the 'phylogenetic.endemism’ function, tests whether observed phylogenetic diversity/endemism is higher than expected, using non-parametric methods
#
#Usage --
#
#For example:
#
#mite.PE <- phylogenetic.endemism(mite, records="site", site.coords=mite.xy, sep.comm.spp="none", phylo.tree=mite.tree, sep.phylo.spp="none", weight.type="geo")
#source(“pe.null.test.R”)
#pe.mite.test <- pe.null.test(mite.PE)
#
#Arguments
#
#phylogenetic.endemism.output : the returned object from a phylogenetic.endemism.R run, a list of length 3 (for PD, pe.type="unweighted") or 7 (for PE: pe.type="weighted"): $PD (numeric vector of phylogenetic diversity/endemism scores for grid cells); $PD_raster (raster layer of PD/PE scores); $grid.matrix (species against cells occurrence matrix - row names match $PD names); $ranges (numeric vector of range weights for phylogenetic branches);  $phyloMatrix (branches against species matrix); $phylo.cell.matrix (branches against cells matrix); $edge.lengths (numeric vector of edge.length in phylo.tree).
#
#nrep : desired number of replicates when generating a null distribution from a random draw of species. Default is 100 for speed (slow for large datasets), but at least 1000 is recommended to ensure smooth null distributions and useful p-values
#
#pe.type="weighted": refers back to how phylogenetic.endemism.R was run; alternatively set to pe.type="unweighted" for phylogenetic diversity with edge length uweighted by range size if this was the original setting (in which case 'phylo.tree' must be supplied).
#
#phylo.tree: the phylogenetic tree that was used to run phylogenetic.endemism.R, only required for testing unweighted pd (pe.type="unweighted").
#
#Details --
#
#With the outputs from the 'phylogenetic.endemism' function, performs the following tests:
#1) non-parametric significance test as to whether observed phylogenetic diversity/endemism is higher or lower than expected, given species richness (and observed species frequencies)
#2) identifies and maps outliers (i.e. in terms of map grid cells that have higher or lower PD/PE) based on quantiles. As categorical: whether score lies more than 1.5 (or other user-defined amount) times outside the interquartile range; as continuous: the factor of the interquartile by which observed values differ from the median / 50% quantile). Returns vectors of values plus raster maps.
#
#Raw phylogenetic diversity/endemism scores are biased both by the completeness of species sampling and species richness itself. This function does not correct the scores per se, but compares them to a null distribution. This is achieved by making replicate random draws from the species pool based on the observed species richness (i.e. same number of species) and the actual species frequencies (more frequent species more likely to be drawn). The distribution of the resulting set of null scores is compared to the observed scores and subsequently grid cells can be mapped as higher or lower than expected (based on significance testing and comparison to null quantiles).
#
#Value --
#
#Returns a list with following elements. Plots of observed and expected phylogenetic diversity/endemism against species richness are generated, as are plots of the generated rasters.
#
#$Quantile.25 : vector of expected lower interquartile range for a given species richness
#
#$Quantile.75 : vector of expected upper interquartile range for a given species richness
#
#$out.above.below : numeric vector assigning outliers categorically, with -1 for lower outlier, 0 non-outlier, 1 for upper outlier, for each grid cell
#
#$out.above.below.raster : map of $out.above.below
#
#$out.continuous : scores for each grid cell = by what factor of the interquartile range observed phylogenetic diversity/endemism differs from the median / 50% quantile
#
#$out.continuous.raster : map of $out.continuous
#
#$P.above : vector of p-values for grid cells having higher than expected phylogenetic diversity/endemism. i.e. very low score = likely to be higher than expected; very high scores = likely to be lower than expected
#
#$P.above.raster : map of $P.above with default 2-tailed colour scheme
#
#$richness : simple observed species richness scores for cells
#
#Required packages --
#
#raster, ape
#
#Authors --
#
#Greg R. Guerin
#
#References --
#
#Guerin, G.R. and Lowe, A.J. (submitted) Mapping phylogenetic endemism in R using georeferenced branch extents.
#
#License --
#
#GNU GPL-3
#
#Version --
#
#1.0
#
{
	require(raster)
	require(ape)
	
	if(class(phylogenetic.endemism.output) != "list") {
		stop("Input data must be a list")
	}
	
	if(pe.type == "weighted") {
		if(!(all(c("PD", "ranges", "PD_raster", "grid.matrix", "phylo.cell.matrix", "edge.lengths", "phyloMatrix") %in% names(phylogenetic.endemism.output)))) {
			stop("Input data list must include elements $PD, $PD_raster, $ranges, $grid.matrix, $phylo.cell.matrix, $edge.lengths, $phyloMatrix; i.e. as output from function phylogenetic.endemism.R")
		} #cls if not all names
		if(!(class(phylogenetic.endemism.output$PD) == "numeric" & class(phylogenetic.endemism.output$PD_raster) == "RasterLayer" & class(phylogenetic.endemism.output$ranges) == "numeric" & class(phylogenetic.endemism.output$grid.matrix) == "data.frame" & class(phylogenetic.endemism.output$phylo.cell.matrix) == "matrix" & class(phylogenetic.endemism.output$edge.lengths) == "numeric" & class(phylogenetic.endemism.output$phyloMatrix) == "matrix")) {
		stop("Input data are not in required format: $PD (numeric), $PD_raster (RasterLayer), $ranges (numeric), $grid.matrix (data.frame), $phylo.cell.matrix (matrix), $edge.lengths (numeric), $phyloMatrix (matrix).")
	} #cls if class not right
	} #cls if weighted
	
	
	if(pe.type == "unweighted") {
		if(!(all(c("PD", "PD_raster", "grid.matrix") %in% names(phylogenetic.endemism.output)))) {
			stop("Input data list must include elements $PD, $PD_raster, $grid.matrix; i.e. as output from function phylogenetic.endemism.R")
		} #cls if not all names
		if(!(class(phylogenetic.endemism.output$PD) == "numeric" & class(phylogenetic.endemism.output$PD_raster) == "RasterLayer" & class(phylogenetic.endemism.output$grid.matrix) == "data.frame")) {
		stop("Input data are not in required format: $PD (numeric), $PD_raster (RasterLayer), $grid.matrix (data.frame)")
	}#cls if not correect object classes
	if(missing(phylo.tree)) {stop("Please provide 'phylo.tree'")}
	} #cls if unweighted
	
				
	cat("Calculating species random sample probabilities and species richness", "\n")
	probs <- colSums(phylogenetic.endemism.output$grid.matrix)/nrow(phylogenetic.endemism.output$grid.matrix) 
	richness <- rowSums(phylogenetic.endemism.output$grid.matrix)[which(rownames(phylogenetic.endemism.output$grid.matrix) %in% names(phylogenetic.endemism.output$PD))]
	
	plot(phylogenetic.endemism.output$PD ~ richness[names(phylogenetic.endemism.output$PD)], pch=20, cex=0.6, main = "PE v species richness (black = observed, red = null)") 
	
	cat("Generating null distributions for range of observed species richness:", "\n")
	n <- 0
	rand.endemism <- list()
	Quantile.75 <- list()
	Quantile.25 <- list()
	Quantile.50 <- list()
	for(i in sort(unique(richness))) {
		n <- n + 1
		rand.endemism[[n]] <- rep(NA, nrep)
		for(j in 1:nrep) {
			
			present <- sample(colnames(phylogenetic.endemism.output$grid.matrix), size = i, prob=probs, replace=FALSE)
			
			
			if(i == 1) {
				rand.endemism[[n]][j] <- 0
			} else {
				
				if(pe.type == "unweighted") {
					phyloDfunc <- function(x,...) {
						x <- as.character(unique(na.omit(x)))
						if(length(x) < 2) {value <- NA} else {
							sub.tree <- drop.tip(phylo.tree, which(!phylo.tree$tip.label %in% x))
							value <- sum(sub.tree$edge.length)
							} #cls else just above...
							return(value)
						} #cls phyloDfunc
					rand.endemism[[n]][j] <- phyloDfunc(present)
					points(i, rand.endemism[[n]][j], col="red", pch=20, cex=0.3)
				} #cls if(pe.type="unweighted
				
				if(pe.type == "weighted") {
					phyloMatrix2 <- phylogenetic.endemism.output$phyloMatrix[which(row.names(phylogenetic.endemism.output$phyloMatrix) %in% present),]
					phyloMatrix2 <- phyloMatrix2[,which(colSums(phyloMatrix2) > 0)]
					phyloMatrix2 <- phyloMatrix2[1,]
					phyloMatrix2[] <- 1
					INVphyloMatrix <- phyloMatrix2 / phylogenetic.endemism.output$ranges[names(phyloMatrix2)]		
					rand.endemism[[n]][j] <- sum(INVphyloMatrix * phylogenetic.endemism.output$edge.lengths[names(phyloMatrix2)])
					points(i, rand.endemism[[n]][j], col="red", pch=20, cex=0.3)
				} #cls if(pe.type="weighted"...
			} #cls else rom high above...
		} #cls for j in nrep
		Quantile.75[[n]] <- quantile(rand.endemism[[n]], na.rm = TRUE)[4]
		Quantile.25[[n]] <- quantile(rand.endemism[[n]], na.rm = TRUE)[2]
		Quantile.50[[n]] <- quantile(rand.endemism[[n]], na.rm = TRUE)[3]
	} #cls for i in sort...
	names(rand.endemism) <- sort(unique(richness))
	names(Quantile.75) <- sort(unique(richness))
	names(Quantile.25) <- sort(unique(richness))
	names(Quantile.50) <- sort(unique(richness))
	
	dev.new()
	plot(phylogenetic.endemism.output$PD ~ richness[names(phylogenetic.endemism.output$PD)], pch=20, cex=0.1, main = "PE v species richness (black = observed, red = null IQR)")
	points(sort(unique(richness)), Quantile.75, type="l", lwd=2, col="red")
	points(sort(unique(richness)), Quantile.25, type="l", lwd=2, col="red")

		
	

	cat("Calculating outliers", "\n")
	out.above.below <- phylogenetic.endemism.output$PD 
	out.continuous <- phylogenetic.endemism.output$PD 

	for (i in names(phylogenetic.endemism.output$PD)) { 
		obs.end <- phylogenetic.endemism.output$PD[which(names(phylogenetic.endemism.output$PD) == i)] 
		assoc.rich <- richness[which(names(richness) == i)] 
		
		upper <- as.numeric(Quantile.75[which(names(Quantile.75) == assoc.rich)])
		lower <- as.numeric(Quantile.25[which(names(Quantile.25) == assoc.rich)])
		mid <- as.numeric(Quantile.50[which(names(Quantile.50) == assoc.rich)])
		inter <- upper - lower 
		up.out <- upper + outlier.range*inter 
		low.out <- lower - outlier.range*inter
		
		if(obs.end > up.out) {whether.out <- 1}
		if(obs.end < low.out) {whether.out <- -1}
		if(obs.end <= up.out & obs.end >= low.out) {whether.out <- 0}
		out.above.below[which(names(out.above.below) == i)] <- whether.out
		
		out.factor.continuous <- (obs.end - mid)/inter 
		out.continuous[which(names(out.continuous) == i)] <- out.factor.continuous
	}
	
	cat("Mapping outlier results", "\n")
	out.above.below.raster <- phylogenetic.endemism.output$PD_raster 
	out.above.below.raster[] <- NA 
	out.above.below.raster[as.numeric(names(out.above.below))] <- out.above.below 
	
	out.breaks <- c(-1.5, -0.5, 0.5, 1.5)
	dev.new()
	plot(out.above.below.raster, breaks=out.breaks, col=heat.colors(3), main = "Outliers - categorical") 
	
	out.continuous.raster <- phylogenetic.endemism.output$PD_raster 
	out.continuous.raster[] <- NA 
	out.continuous.raster[as.numeric(names(out.continuous))] <- out.continuous 

	dev.new()
	plot(out.continuous.raster, main = "Outliers - continuous")


	cat("Calculating significance", "\n")
	P.above <- phylogenetic.endemism.output$PD
	for (i in names(phylogenetic.endemism.output$PD)) {
		obs.end <- phylogenetic.endemism.output$PD[which(names(phylogenetic.endemism.output$PD) == i)]
		assoc.rich <- richness[which(names(richness) == i)]
		randoms <- rand.endemism[which(names(rand.endemism) == assoc.rich)]
		randoms <- sapply(randoms, FUN=function(x) {as.numeric(x)})[,1]
		P.above[which(names(P.above) == i)] <- (1 + length(randoms[randoms >= obs.end])) / (1 + length(randoms))
	}
	
	cat("Mapping significance results", "\n")
	P.above.raster <- phylogenetic.endemism.output$PD_raster 
	P.above.raster[] <- NA 
	P.above.raster[as.numeric(names(P.above))] <- P.above

	p.breaks <- c(0, 0.0001, 0.001, 0.01, 0.025, 0.975, 0.99, 0.999, 0.9999, 1)
	dev.new()
	plot(P.above.raster, breaks=p.breaks, col=topo.colors(9), main = "Outliers - significance (e.g. higher <0.025; lower >0.975")
	
	
	cat("Collating outputs", "\n")
	outputs <- list(Quantile.25 = Quantile.25, Quantile.75 = Quantile.75, out.above.below = out.above.below, out.above.below.raster = out.above.below.raster, out.continuous = out.continuous, out.continuous.raster = out.continuous.raster, P.above = P.above, P.above.raster = P.above.raster, richness = richness[names(phylogenetic.endemism.output$PD)])
	return(outputs)

			
			
}#cls function


