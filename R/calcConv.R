#'Quantifies convergent evolution by the ratio of the current to maximum past phenotypic distance between two or 
#'more lineages, as described in Stayton (2015). Code written by Jonathan S. Mitchell for Zelditch et al. (2017).
#'
#'calcConv prepares arguements for the CalcCs function
#'
#'@param phy The phylogeny of interest in phylo format
#'@param traits Phenotypic data for all tips
#'@param focaltaxa A list consisting of the names of all putatively convergent taxa
#'@param anc A matrix of user supplied ancestral trait values at internal nodes (formatted as "traits" but with node number as rownames)
#'@param VERBOSE Whether or not to print progress
#'
#'@details calcConv is a wrapper function which formats data, performs ancestral state reconstructions, obtains 
#'distance matrices, and determines pairwise combinations of focal taxa, which are then used as arguements for
#'the CalcCs function, which calculates values for C1-C4 of Stayton (2015) for each pair of putatively 
#'convergent taxa. 
#'
#'@return C1-C4 convergence measures for all pairs of putatively convergent taxa.
#'
#'@import ape phytools
#'
#'@importFrom stats dist
#'@importFrom utils combn
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
#'answer<-calcConv(phy,traits,focaltaxa,anc=NULL,VERBOSE=FALSE)

calcConv <- function (phy, traits, focaltaxa, anc = NULL, VERBOSE = FALSE) 
{
    if (!inherits(phy, "phylo")) {
        stop("first argument, phy, must be of class phylo!")
    }
    if (ncol(as.matrix(traits)) == 1) {
        traits <- cbind(traits, 0)
    }
    if (!inherits(traits, "matrix")) {
        Names <- names(traits)
        traits <- matrix(traits, ncol = 1)
        rownames(traits) <- Names
    }
    if (nrow(traits) != Ntip(phy)) {
        stop("Number of rows does not match number of tips!")
    }
    if (length(intersect(rownames(traits), phy$tip.label)) != 
        Ntip(phy)) {
        print("Warning: the row names in traits do not match tip labels in phy; assume that they are in the same order.")
        rownames(traits) <- phy$tip.label
    }
    traits <- traits[phy$tip.label, ]
    if (VERBOSE == TRUE) {
        print("reached traits")
    }
    
	if(is.null(anc)){								## Changed block below to pass user defined ancestral states if object anc is not NULL
    ACE <- apply(traits, 2, fastAnc, tree = phy)	## or if anc is not provided by the user calculate ancestral states as in previous versions
	} else {										##
	ACE <- anc										##
	print("using user supplied ancestral states")	##
	}												##

    tipNums <- as.character(1:Ntip(phy))
    names(tipNums) <- phy$tip.label
    focaltaxa <- tipNums[focaltaxa]
    rownames(traits) <- tipNums[rownames(traits)]
    allVals <- rbind(traits, ACE)
    allDists <- as.matrix(dist(allVals, method = "euclidean"))
    Ancs <- lapply(tipNums, pullNodeSeq, phy = phy)
    names(Ancs) <- tipNums
    Combinations <- combn(focaltaxa, 2)
    if (VERBOSE == TRUE) {
        print("starting combinations...")
    }
    Cmat <- apply(Combinations, 2, function(x) calcCs(tips = x[1:2], 
        ancList = Ancs, allDists = allDists, phy = phy, VERBOSE = VERBOSE))
    Combinations <- Combinations[1:2, ]
    out <- apply(Cmat, 1, mean)
    return(out)
}