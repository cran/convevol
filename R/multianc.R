#'Reconstructs ancestral states for multiple characters
#'
#'Uses fastAnc to reconstruct ancestral states for multiple phenotypic characters
#'
#'@param phy The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'
#'@details None
#'
#'@return A matrix with the tips data in the first n rows and the ancestral data in the remaining n-1 rows.  
#'
#'@import ape MASS phytools
#'
#'@export
#'
#'@references Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics
#'and evolution in R langauge. Bioinformatics, 20, 289-290.
#'
#'Revell, L. J. (2012) phytools: An R package for phylogenetic comparative 
#'biology (and other things). Methods Ecol. Evol. 3 217-223.
#'
#'@examples
#'
#'phy<-rtree(10)
#'
#'phendata<-fastBM(phy,nsim=2)
#'
#'ancs<-multianc(phy,phendata)



multianc<-function(phy,phendata)

{

#Error checking

if (!inherits(phy,"phylo")) 
	stop("your tree must be class 'phylo.'")

if (nrow(phendata) != length(phy$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    
if (is.null(rownames(phendata))) {
	warning("no row names for data.  Assuming that the rows are in the same order as tips.")
      rownames(X) <- phy$tip.label
	}

#The function

firstvar<-fastAnc(phy,phendata[, 1])

allancstates<-matrix(data=0, length(firstvar), ncol(phendata))

for (i in 1:ncol(phendata)) {allancstates[, i]<-fastAnc(phy,phendata[, i])}

colnames(allancstates)<-colnames(phendata)

alldata<-rbind(phendata,allancstates)

}

