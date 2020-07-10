phylogenetic.diversity.sites <- function(species_records, phylo.tree)
{


	if(class(species_records) != "data.frame") {
		stop("Species data must be in a data.frame")
	}


	if(class(phylo.tree) != "phylo") {
		stop("Phylogenetic tree must be in 'phylo' format")
	}


		cat("Checking that occurrence matrix and phylo match." , "\n")
		#trim excess species in species_records that aren't in phylo.tree
		species_records <- species_records[,which(colnames(species_records) %in% phylo.tree$tip.label)]
		#trim excess tips in phylo.tree that aren't in species_records
		phylo.tree <- drop.tip(phylo.tree, which(!(phylo.tree$tip.label %in% colnames(species_records))))



		cat("Calculating phylogenetic diversity", "\n")

		PDmatFunc <- function(x, phylo.tree) {
			x <- x[,which(x>0)]
			drops <- setdiff(phylo.tree$tip.label, colnames(x))
			if(class(x) == "data.frame") {sub <- ape::drop.tip(phylo.tree, drops)
			return(sum(sub$edge.length))} else {return(NA)}
			}

		PD <- list()
		for(i in 1:nrow(species_records)) {
			cell_row <- species_records[i,]
			PD[i] <- PDmatFunc(x=cell_row, phylo.tree=phylo.tree)
			}
		PD <- unlist(PD)
		names(PD) <- row.names(species_records)
		PD <- PD[!is.na(PD)]


		cat("Collating outputs", "\n")
		output <- list()
		output$PD <- PD

		return(output)


	} #cls function
