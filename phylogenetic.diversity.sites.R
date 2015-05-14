phylogenetic.diversity.sites <- function(species_records, phylo.tree)
#
#Description --
#
#Calculates (unweighted and uncorrected) phylogenetic diversity for a set of sample sites.
#
#Usage --
#
#For example:
#####Preparation for this example:#library(vegan) #load the vegan package#data(mite)#library(ape) #for rtree function:#mite.tree <- rtree(n=ncol(mite), tip.label=colnames(mite)) #for this example, generate a phylogenetic tree of the species in the mite dataset with random relationships and branch lengths####Usage of the function:
#source(“phylogenetic.diversity.sites.R”)#mite.PD <- phylogenetic.diversity.sites(mite, phylo.tree=mite.tree)
#
#Arguments --
#
#species_records: a data.frame with rows as sites and columns as species
#
#phylo.tree: a phylogenetic tree of class 'phylo' containing species in the occurrence data. Must have branch lengths.
#
#Details --
#
#Unweighted and uncorrected phylogenetic diversity for community samples organised into sites.
#
#Value --
#
#Returns a vector of PD scores (list with $PD as the vector) for sites in the community matrix
#
#
#Required packages --
#
#ape
#
#Authors --
#
#Greg R. Guerin
#
#References --
#
#[code modified from:–] Guerin, G.R. and Lowe, A.J. (submitted) Mapping phylogenetic endemism in R using georeferenced branch extents.
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
		require(ape)
	
	
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
			if(class(x) == "data.frame") {sub <- drop.tip(phylo.tree, drops)
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