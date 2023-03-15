#'Extracts a vector of ancestors for a given taxon. Code written by Jonathan S. Mitchell for 
#'Zelditch et al. (2017)
#'
#'pullNodeSeq Extracts a vector of all ancestors of a given taxon in a phylogeny.    
#'
#'@param phy The phylogeny of interest in phylo format
#'@param tip The tip of interest
#'
#'@return A vector of ancestors  
#'
#'@import MASS phytools
#'
#'@export
#'
#'@references Zelditch, M.L., J. Ye, J.S. Mitchell, and D.L. Swiderski. 2017. Rare ecomorphological
#'convergence on a complex adaptive landscape: Body size and diet mediate evolution of jaw shape in
#'squirrels (Sciuridae). Evolution 71:633-649.
#'
#'@examples
#'
#'phylogeny<-rtree(100)
#'answer<-pullNodeSeq(phy=phylogeny,tip="t1")

pullNodeSeq <- function(phy, tip)	{
	if (is.na(as.numeric(tip)))	{
		tip <- which(phy$tip.label == tip)
	}
	tip <- as.numeric(tip)
	Stop <- Ntip(phy) + 1
	Anc <- phy$edge[which(phy$edge[,2] == tip),1]
	if (Anc != Stop)		{
		Anc <- c(Anc, pullNodeSeq(phy, Anc))
	}
	return(as.character(Anc))
}