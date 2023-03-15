#'Uses simulations to assess the significance of C1-C4 measures of convergent 
#'evolution as described in Stayton (2015). Code written by Jonathan S. 
#'Mitchell for Zelditch et al. (2017)
#'
#'convSig  calculates the significance of measures of convergent evolution.
#'
#'@param phy The phylogeny of interest in phylo format
#'@param traits Phenotypic data for all tips
#'@param focaltaxa A list consisting of the names of all putatively convergent taxa
#'@param nsim The number of simulations to use to assess significance
#'
#'@details This script simulates data according to a Brownian motion model of
#'evolution, and then assesses convergene on that simulated data. The number 
#'of times that the simulated data produces greater convergence than that 
#'observed in the empirical data is used to calculate a P-value. 
#'
#'@return C1-C4 convergence measures for all pairs of putatively convergent 
#'taxa and their associated p-values.
#'
#'@import ape MASS phytools
#'
#'@importFrom geiger sim.char
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
#'@examples
#'
#'phy<-rtree(100)
#'traits<-fastBM(phy,nsim=3)
#'focaltaxa<-c("t1","t50","t100")
#'answer<-convSig(phy,traits,focaltaxa,nsim=10)

convSig <- function(phy, traits, focaltaxa, nsim=1e3)	{
	data <- calcConv(phy, traits, focaltaxa)
	
	# run simulations
	phylMat <- vcv.phylo(phy)
	phylMat2 <- phyl.vcv(traits, phylMat, 0)
	simDat <- sim.char(phy, phylMat2$R, nsim, model="BM", root=0)
	simOut <- apply(simDat, 3, calcConv, phy=phy, focaltaxa=focaltaxa)
	pvals <- sapply(1:4, function(x) length(which(simOut[x,] >= data[x]))) / nsim
	out <- cbind(data, pvals)
	return(out)
}