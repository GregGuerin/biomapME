endemism.null.test <- function(weighted.endemism.output, nrep = 100, outlier.range = 1.5)
#
#Description --
#
#Taking the outputs from the 'weighted.endemism' function, tests whether observed endemism is higher than expected, using non-parametric methods
#
#Usage --
#
#For example:
#
#endemism_mydata <- weighted.endemism(mite, site.coords=mite.xy, records="site") 
#
#endemism.test.example <- endemism.null.test(endemism_mydata)
#
#Arguments
#
#weighted.endemism.output : the returned object from a weighted.endemism run (or equivalent constructed separately). Consists of a list of length 4: $WE (named vector of endemism scores for grid cells), $WE_raster (raster layer of endemism scores), $weights (named vector of range weights for species used in the calculation of weighted endemism); $grid.matrix (species against cells occurrence matrix - row names must match $WE names, column names must match $weights names). If 'corrected weighted endemsim' was specified in the weighted.endemism run, this function will not accept the output. This is because correcting via dividing endemism by species richness is an alternative to this function, which uses more sophisticated methods for testing whether endemism is different than that expected for a given species richness.
#
#nrep : desired number of replicates when generating a null distribution from a random draw of species. Default is 100 for speed, but at least 1000 is recommended to ensure smooth null distributions and useful p-values
#
#Details --
#
#With the outputs from the 'weighted.endemism' function, performs the following tests:
#1) non-parametric significance test as to whether observed endemism is higher or lower than expected, given species richness (and observed species frequencies)
#2) identifies and maps outliers (i.e. in terms of map grid cells that have higher or lower endemism) based on quantiles. As categorical: whether endemism score lies more than 1.5 (or other user-defined amount) times outside the interquartile range; as continuous: the factor of the interquartile by which observed values differ from the median / 50% quantile). Returns vectors of values plus raster maps.
#
#Raw weighted endemism scores are biased both by the completeness of species sampling and species richness itself. Correcting by dividing by the observed number of species ('corrected weighted endemism' of Crips et al. 2001) is a proposed correction, but the relationship between endemism scores and species richness is not linear under a null model (random species draws), as increasingly infrequent species are drawn as richness increases, thereby increasing CWE. While correcting endemism scores in a more sophisticaed way is possible, this function does not correct the scores per se, but compares them to a null distribution. This is achieved by making replicate random draws from the species pool based on the observed species richness (i.e. same number of species) and the actual species frequencies (more frequent species more likely to be drawn). The distribution of the resulting set of null endemsim scores is compared to observed endemism and subsequently grid cells can be mapped as higher or lower than expected (based on significance testing and comparison to null quantiles).
#
#Value --
#
#Returns a list with following elements. Plots of observed and expected endemism against species richness are generated, as are plots of the generated rasters.
#
#$Quantile.25 : vector of expected lower interquartile range for a given species richness
#
#$Quantile.75 : vector of expected upper interquartile range for a given species richness
#
#$out.above.below : vector assigning outliers categorically, with -1 for lower outlier, 0 non-outlier, 1 for upper outlier for each grid cell
#
#$out.above.below.raster : map of $out.above.below
#
#$out.continuous : scores for each grid cell = by which factor of the interquartile range observed endemism differs from the median / 50% quantile
#
#$out.continuous.raster : map of $out.continuous
#
#$P.above : vector of p-values for grid cells having higher than expected endemism. Very low scores = more likely to be higher than expected; very high scores = more likely to be lower than expected
#
#$P.above.raster : map of $P.above
#
#$richness : simple observed species richness scores for cells
#
#Required packages --
#
#raster
#
#Authors --
#
#Greg R. Guerin & Lasse Ruokolainen
#
#References --
#
#Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12361
#
#License --
#
#GPL-3
#
{
	
	require(raster)
	
	if(class(weighted.endemism.output) != "list") {
		stop("Input data must be a list")
	}
	if(!(all(c("WE", "WE_raster", "weights", "grid.matrix") %in% names(weighted.endemism.output)))) {
		stop("Input data list must include elements $WE, $WE_raster, $weights, $grid.matrix, i.e. as output from function weighted.endemism")
	}
	if(!(class(weighted.endemism.output$WE) == "numeric" & class(weighted.endemism.output$WE_raster) == "RasterLayer" & class(weighted.endemism.output$weights) == "numeric" & class(weighted.endemism.output$grid.matrix) == "data.frame")) {
		stop("Input data are not in required format: $WE (numeric), $WE_raster (RasterLayer), $weights (numeric), $grid.matrix (data.frame)")
	}
				
	cat("Calculating species random sample probabilities and species richness", "\n")
	probs <- colSums(weighted.endemism.output$grid.matrix)/nrow(weighted.endemism.output$grid.matrix) 
	richness <- rowSums(weighted.endemism.output$grid.matrix) #vector of observed spp richness for each grid cell
	
	plot(weighted.endemism.output$WE ~ richness, pch=20, cex=0.6, main = "Endemism against species richness (black = observed, red = null)") 
	
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
			present <- sample(colnames(weighted.endemism.output$grid.matrix), size = i, prob=probs, replace=FALSE) 
			ranges_subset <- weighted.endemism.output$weights[which(names(weighted.endemism.output$weights) %in% present)]
			rand.endemism[[n]][j] <- sum(1  / ranges_subset)
			points(i, rand.endemism[[n]][j], col="red", pch=20, cex=0.3)
		}
		Quantile.75[[n]] <- quantile(rand.endemism[[n]])[4]
		Quantile.25[[n]] <- quantile(rand.endemism[[n]])[2]
		Quantile.50[[n]] <- quantile(rand.endemism[[n]])[3]
	}
	names(rand.endemism) <- sort(unique(richness))
	names(Quantile.75) <- sort(unique(richness))
	names(Quantile.25) <- sort(unique(richness))
	names(Quantile.50) <- sort(unique(richness))
	
	dev.new()
	plot(weighted.endemism.output$WE ~ richness, pch=20, cex=0.1, main = "Endemism against species richness (black = observed, red = null IQR)")
	points(sort(unique(richness)), Quantile.75, type="l", lwd=2, col="red")
	points(sort(unique(richness)), Quantile.25, type="l", lwd=2, col="red")
		
	
	#####OUTLIERS - categorical and continuous
	cat("Calculating outliers", "\n")
	out.above.below <- weighted.endemism.output$WE 
	out.continuous <- weighted.endemism.output$WE 

	for (i in names(weighted.endemism.output$WE)) { 
		obs.end <- weighted.endemism.output$WE[which(names(weighted.endemism.output$WE) == i)] 
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
	out.above.below.raster <- weighted.endemism.output$WE_raster 
	out.above.below.raster[] <- NA 
	out.above.below.raster[as.numeric(names(out.above.below))] <- out.above.below 
	
	out.breaks <- c(-1.5, -0.5, 0.5, 1.5)
	dev.new()
	plot(out.above.below.raster, breaks=out.breaks, col=heat.colors(3), main = "Outliers - categorical") 
	
	out.continuous.raster <- weighted.endemism.output$WE_raster 
	out.continuous.raster[] <- NA 
	out.continuous.raster[as.numeric(names(out.continuous))] <- out.continuous 

	dev.new()
	plot(out.continuous.raster, main = "Outliers - continuous")


	####SIGNIFICANCE
	cat("Calculating significance", "\n")
	P.above <- weighted.endemism.output$WE 
	for (i in names(weighted.endemism.output$WE)) {
		obs.end <- weighted.endemism.output$WE[which(names(weighted.endemism.output$WE) == i)]
		assoc.rich <- richness[which(names(richness) == i)]
		randoms <- rand.endemism[which(names(rand.endemism) == assoc.rich)]
		randoms <- sapply(randoms, FUN=function(x) {as.numeric(x)})[,1]
		P.above[which(names(P.above) == i)] <- (1 + length(randoms[randoms >= obs.end])) / (1 + length(randoms))
	}
	
	cat("Mapping significance results", "\n")
	P.above.raster <- weighted.endemism.output$WE_raster 
	P.above.raster[] <- NA 
	P.above.raster[as.numeric(names(P.above))] <- P.above

	p.breaks <- c(0, 0.0001, 0.001, 0.01, 0.025, 0.975, 0.99, 0.999, 0.9999, 1)
	dev.new()
	plot(P.above.raster, breaks=p.breaks, col=topo.colors(9), main = "Outliers - significance (e.g. higher than expected <0.025; lower >0.975")
	
	
	cat("Collating outputs", "\n")
	outputs <- list(Quantile.25 = Quantile.25, Quantile.75 = Quantile.75, out.above.below = out.above.below, out.above.below.raster = out.above.below.raster, out.continuous = out.continuous, out.continuous.raster = out.continuous.raster, P.above = P.above, P.above.raster = P.above.raster, richness = richness)
	return(outputs)
	
}


