#'Calculates the C1-C4 measures of convergent evolution between two lineages as described in 
#'Stayton (2015). All measures quantify convergence by the ratio of current to maximum past phenotypic
#'distance between lineages. Can be used as-is but more often will be used within the calcConv script.
#'Code written by Jonathan S. Mitchell for Zelditch et al. (2017)
#'
#'calcCs calculates the C1-C4 measures of convergent evolution
#'
#'@param tips Two putatively convergent tips
#'@param ancList A list of ancestors of all tips. Most often obtained from calcConv 
#'@param allDists A matrix of phenotypic distances between all nodes (tips and ancestors). Most often
#'obtained from calcConv
#'@param phy The phylogeny of interest
#'@param VERBOSE Whether or not to print progress
#'
#'@details calcCs calculates values of C1-C4, all of which are fundamentally based on comparing the 
#'current phenotypic distance between two tips to the maximum past distances between the ancestors of 
#'those tips. Higher values indicate a greater amount of past phenotypic distance which has been "closed"
#'by subsequent evolution, and thus greater convergence. C1 is the ratio of tip to maximum ancestral 
#'distance. C2 is the difference of those two values. C3 scales C2 by the total amount of evolution that
#'has occured in the two lineages. C4 scales C2 by the total amount of evolution that has occurred in the
#'entire phylogeny. The arguments for this function will usually be obtained from the calcConv script in
#'convevol - this allows certain computationally-intensive steps (e.g., calculating ancestral states) to
#'only be performed once, thus saving a great deal of time. This script also corrects an error in the 
#'calculation of C4 of previous versions of convevol.   
#'
#'@return C1-C4 convergence measures for all pairs of putatively convergent taxa.
#'
#'@import phytools
#'
#'@export
#'
#'@references Stayton, C.T. 2015. The definition, recognition, and interpretation of convergent evolution, and 
#'two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140-2453.
#'
#'Zelditch, M.L., J. Ye, J.S. Mitchell, and D.L. Swiderski. 2017. Rare ecomorphological
#'convergence on a complex adaptive landscape: Body size and diet mediate evolution of jaw shape in
#'squirrels (Sciuridae). Evolution 71:633-649.
#'
#'
calcCs <- function(tips, ancList, allDists, phy, VERBOSE=FALSE)	{
	if (VERBOSE == TRUE)		{
		cat("Analzying ", tips)
	}
	# tip distances
	Dtip <- allDists[tips[1], tips[2]]
	
	# what nodes are shared between tip1 and tip2?
	overLap <- intersect(ancList[[tips[1]]], ancList[[tips[2]]])
	
	# what node is the mrca of tip1 and tip2?
	MRCA <- max(as.numeric(overLap))
	
	# what nodes are UNIque to tip 1 and what are UNIque to tip 2? E.g., what nodes lead from the MRCA to each tip?
	unit1 <- c(setdiff(c(tips[1], ancList[[tips[1]]]), overLap), MRCA)
	unit2 <- c(setdiff(c(tips[2], ancList[[tips[2]]]), overLap), MRCA)
	
	# compare all post-MRCA nodes along each lineage to one another (including the tips)
	ancPairs <- sapply(unit1, function(t1) sapply(unit2, function(t2) allDists[t1,t2]))
	
	# sum changes along each lineage including the MRCA
	lin1 <- c(unit1, as.character(MRCA))
	Lin1 <- sum(sapply(1:(length(lin1)-1), function(x) allDists[lin1[x],lin1[x+1]]))
	lin2 <- c(unit2, as.character(MRCA))
	Lin2 <- sum(sapply(1:(length(lin2)-1), function(x) allDists[lin2[x],lin2[x+1]]))
	
	# sum all pairwise distances between nodes of the subtree from the MRCA of tip1 & tip2
	totalClade <- as.character(getDescendants(phy, MRCA))
	totalMove <- sum(allDists[totalClade, totalClade]) / 2
	
	# Dmax as the maximum distance between ancestor pairs
	Dmax <- max(ancPairs)
	C1 <- 1 - (Dtip / Dmax)
	C2 <- Dmax - Dtip
	C3 <- C2 / (Lin1 + Lin2)
	C4 <- C2 / totalMove
	return(c(C1=C1, C2=C2, C3=C3, C4=C4))
}