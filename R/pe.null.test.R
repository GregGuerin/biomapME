pe.null.test <- function(phylogenetic.endemism.output, nrep = 100, outlier.range = 1.5, pe.type="weighted", phylo.tree)
  {

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


