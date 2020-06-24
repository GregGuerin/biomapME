endemism.null.test <- function(weighted.endemism.output, nrep = 100, outlier.range = 1.5)

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


