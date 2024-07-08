#'Calculates the number of measurements that can be made between two lineages for each pairwise comparison within a set of putatively convergent tips (group identity may also be taken into account). Useful for determining which comparisons are not informative, and constructing a group object before running calcConvCt or convSigCt.
#'
#'pwCheck Calculates the number of measurements that can be made between two lineages for each pairwise comparison within a set of putatively convergent tips (group identity may also be taken into account). Useful for determining which comparisons are not informative, and constructing a group object before running calcConvCt or convSigCt.
#'
#'@param phy The time calibrated phylogeny of interest in phylo format
#'@param focaltaxa a vector of tip labels for the putatively convergent taxa to be compared
#'@param groups an optional vector of groups with names matching focaltaxa. Indicates the group identity of all putatively convergent taxa and limits Ct measures to intergroup comparisons only
#'@param conservative logical value indicating whether candidate nodes for measurement of Dmax.t should be restricted to occurr before the oldest stem lineage of the two groups involved in each pairwise comparison. Stem lineage age for each group is defined as the height of the parent node of the branch subtending the most recent common ancestor of tips within a group. Where groups include a single tip, the parent node of the tip's subtending branch is used. Requires group object to be provided by user.
#'
#'@return A list of the following components:
#'@return taxa a matrix of uninformative tip comparisons
#'@return path a vector with the number of measurements for all pairwise comparisons - a summary of this is also printed when running the function
#'
#'@import geiger phytools
#'
#'@export
#'
#'@references Grossnickle DM, Brightly WH, Weaver LN, Stanchak KE, Roston RA, Pevsner SK, Stayton CT, Polly PD, Law CJ. 2022. A cautionary note on quantitative measures of phenotypic convergence. in revision
#'Zelditch ML, Ye J, Mitchell JS, Swiderski DL. 2017. Rare ecomorphological convergence on a complex adaptive landscape: Body size and diet mediate evolution of jaw shape in squirrels (Sciuridae). Evolution 71: 633-649
#'Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution and two new measures for quantifying and assessing the significance of convergence. Evolution 69(8): 2140-2153.
#'Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.
#'Felsenstein, J. 1985. Phylogenies and the comparative method. American Naturalist, 125, 1-15.
#'
#'@examples
#'\dontrun{
#'library(phytools)
#'library(geiger)
#'
#'# create time calibrated tree
#'mytree<-rtree(100)
#'mycalibration <- makeChronosCalib(mytree, node="root", age.max=50)
#'phy <- chronos(mytree, lambda = 1, model = "correlated", calibration = mycalibration, 
#'control = chronos.control() )
#'class(phy)<-"phylo"
#'
#'# create three normally distributed phenotypic traits
#'traits <- cbind(rnorm(Ntip(phy)),rnorm(Ntip(phy)),rnorm(Ntip(phy)))
#'rownames(traits) <- phy$tip.label
#'
#'#	select two random tips, excluding sister taxa
#'pairs <- apply(combn(phy$tip.label,2),2,function(x) nodepath(phy,which(phy$tip.label == x[1]),
#'which(phy$tip.label == x[2])))
#'nosis <- combn(phy$tip.label,2)[,unlist(lapply(pairs, function(x) length(x) > 3))]
#'focaltaxa <- nosis[,sample(1:ncol(nosis),1)]
#'
#'pwCheck(phy,focaltaxa)
#'}

pwCheck <- function(phy, focaltaxa, groups = NULL, conservative = FALSE) {
	Combinations <- combn(focaltaxa, 2)
	if(!is.null(groups)){
	Combinations <- Combinations[,apply(Combinations,2,function(x) length(unique(groups[x]))>1)]
	}
	
	lng_path <-
	apply(Combinations,2,function(y) {
	sum(unlist(lapply(which(phy$tip.label %in% y),function(x) length(nodepath(phy,getMRCA(phy,y),x)) - 2)))
	})

if(conservative == TRUE & is.null(groups)){
print("Conservative analyses require a group vector to work correctly, please supply this and rerun")
}
if(conservative == TRUE & !is.null(groups)){
            lim.list <- vector("list", ncol(Combinations))
            for (j in 1:ncol(Combinations)) {
                focal_grps <- groups[names(groups) %in% Combinations[,j]]
                focal_mrca <- lapply(focal_grps, function(x) {
                  mmm <- getMRCA(phy, names(groups[groups == x]))
				  if (is.null(mmm)) {
                   mmm <- max(nodepath(phy, which(phy$tip == names(groups[groups == x])), Ntip(phy) + 1))
                  } else{ mmm <- nodepath(phy,mmm,Ntip(phy) + 1)[2] }
                  return(mmm)
                })
				edge <- as.data.frame(cbind(phy$edge, nodeHeights(phy)))
				colnames(edge) <- c("parent", "child", "parent.height", "child.height")
                focal_height <- min(edge[edge$parent %in% unlist(focal_mrca), ]$parent.height)
                lim.list[[j]] <- focal_height
            }

lng_path <- vector("list",ncol(Combinations))
for(y in 1:ncol(Combinations)){
fcl.path <- lapply(which(phy$tip.label %in% Combinations[,y]),function(x) nodepath(phy,getMRCA(phy,Combinations[,y]),x))
tmp.sum <- sum(unlist(lapply(fcl.path, function(z) 
length(z[unlist(lapply(z,function(x) nodeheight(phy,x))) < lim.list[y]]) - 1
)))
if(tmp.sum < 0){lng_path[[y]] <- 0} else {lng_path[[y]] <- tmp.sum}
}
lng_path <- unlist(lng_path)
}

	names(lng_path) <- apply(Combinations,2,function(x) paste(x,collapse = "-"))
	
	print(summary(lng_path))
	if(!is.null(groups)){
	zero <- Combinations[,which(lng_path == 0)]
	if(is.null(dim(zero))){zero <- matrix(zero)}
	return(list("taxa" = zero,
				"groups" = apply(zero,2,function(x) groups[x]),
				"path" = lng_path
				))
	} else { return(list("taxa" = Combinations[,which(lng_path == 0)],
							"path" = lng_path
							))}
}







