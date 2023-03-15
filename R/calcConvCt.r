#'Computes Ct-metric scores for putatively convergent tips (or groups of tips) given a set of user provided phenotypic characters and a time calibrated phylogeny.
#'
#'calcConvCt Computes Ct-metric scores for putatively convergent tips (or groups of tips) given a set of user provided phenotypic characters and a time calibrated phylogeny.
#'
#'@param phy The time calibrated phylogeny of interest in phylo format
#'@param traits A matrix of numeric phenotypic traits with rownames matching tip labels of phy
#'@param focaltaxa A vector of tip labels for the putatively convergent taxa to be compared
#'@param groups An optional vector of groups with names matching focaltaxa. Indicates the group identity of all putatively convergent taxa and limits Ct measures to intergroup comparisons only
#'@param conservative Logical value indicating whether measurement of Dmax.t should be restricted to before the origin of the oldest lineage in each pairwise comparison of the focaltaxa. The origin of each convergent lineages is taken as the most recent common ancestor of tips in each user defined group. Where groups include a single tip, the parent node of the tip's subtending branch is used. Requires group object to be provided by user.
#'@param VERBOSE Logical value indicating whether model information should be printed during computation
#'
#'@details Function incorporates the optimizations introduced by Zelditch et al. (2017), which significantly improve runtimes
#'@details Reconstructions part way along branches are obtained using equation [2] of Felsenstein (1985), following code modified from the phytools (Revell, 2012) function contMap
#'
#'@return A list of the following components:
#'@return mean A named vector of Ct-metrics averaged from all pairwise comparisons of focaltaxa. If user provided groups, this is based only on comparisons between taxa belonging to different groups.
#'@return Cmat A matrix of Ct-metrics for each pairwise comparison.
#'@return path_df A list of dataframes, one per pairwise comparison of the focal taxa, each containing information from all timepoint measurements of the two putatively convergent lineages. These provide the nodes at which comparisons were drawn, the evolutionary path along which that node fell (i.e., leading to one of two tips), the node height, reconstructed ancestral states at that node for each phenotypic trait, reconstructed ancestral values for each trait along the opposite path, and the phenotypic distance between the two lineages at that point.
#'@return grp.mean A matrix of Ct-metrics summarized for inter-group comparisons, returned only if user defined groups were specified. Provides overall results matching those reported in "mean", results for each unique inter-group comparison, and results averaged with equal weight given to each unique inter-group comparison (i.e., ignoring differences in the number of tips in each group).
#'
#'@import MASS phytools ape
#'
#'@importFrom utils combn
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
#'\donttest{# create time calibrated tree
#'mytree<-rtree(100)
#'mycalibration <- makeChronosCalib(mytree, node="root", age.max=50)
#'phy <- chronos(mytree, calibration = mycalibration, control = chronos.control() )
#'class(phy)<-"phylo"
#'
#'# create three normally distributed phenotypic traits
#'traits <- cbind(rnorm(Ntip(phy)),rnorm(Ntip(phy)),rnorm(Ntip(phy)))
#'rownames(traits) <- phy$tip.label
#'focaltaxa <- sample(phy$tip.label, 5)
#'
#'system.time(run <- calcConvCt(phy, traits, focaltaxa))}

calcConvCt <- function(phy, traits, focaltaxa, groups = NULL, conservative = FALSE, VERBOSE=FALSE)	{
	if (!inherits(phy,"phylo"))	{
		stop("first argument, phy, must be of class phylo!")
	}
	if (ncol(as.matrix(traits)) == 1){traits<-cbind(traits,0)}
	if (!inherits(traits, "matrix"))	{
		Names <- names(traits)
		traits <- matrix(traits, ncol=1)
		rownames(traits) <- Names
	}
	if (nrow(traits) != Ntip(phy))	{
		stop("Number of rows does not match number of tips!")
	}
	if (length(intersect(rownames(traits), phy$tip.label)) != Ntip(phy))		{
		print("Warning: the row names in traits do not match tip labels in phy; assume that they are in the same order.")
		rownames(traits) <- phy$tip.label
	}
	traits <- traits[phy$tip.label,]
	if (VERBOSE == TRUE)		{
		print("reached traits")
	}
	# for C1: need phylogeny, trait matrix, ancestral recos
	ACE <- apply(traits, 2, fastAnc, tree=phy)

	# number the tips
	tipNums <- as.character(1:Ntip(phy))
	names(tipNums) <- phy$tip.label
	focaltaxa <- tipNums[focaltaxa]
	rownames(traits) <- tipNums[rownames(traits)]

	# matrix with node and tip vectors
	allVals <- rbind(traits, ACE)
	
	# edge matrices with ancestral state values for each trait (adapted from contMAP - Revell, 2012)
	A<-vector("list",ncol(allVals))
	for(i in 1:ncol(allVals)){
	A[[i]] <- matrix(allVals[as.character(phy$edge),i], nrow(phy$edge), ncol(phy$edge))
	}
	
	edge<-as.data.frame(cbind(phy$edge,nodeHeights(phy)))
	colnames(edge)<-c("parent","child","parent.height","child.height")
	edge<-lapply(A,function(x) cbind(edge,x))	##rename test
	edge<-lapply(edge,function(x) {colnames(x)[5:6]<-c("parent.anc","child.anc")
									return(x)})



	# get distances between tips and nodes
	allDists <- as.matrix(dist(allVals, method="euclidean"))

	# What node is ancestral to each tip?
	Ancs <- lapply(tipNums, pullNodeSeq, phy=phy)
	names(Ancs) <- tipNums

	# unique combinations of focal taxa, if groups object exists remove within group comparisons
	Combinations <- combn(focaltaxa, 2)
	
	if(!is.null(groups)){
	names(groups)<-tipNums[names(groups)]
	Combinations <- Combinations[,apply(Combinations,2,function(x) length(unique(groups[x]))>1)]
	
	## add row with height of oldest group mrca for each pair of taxa to Combinations matrix (if conservative Dmax option is chosen) - does this fail if only 1 taxon is in one of the focal groups?
	if(conservative == TRUE){
		if (VERBOSE == TRUE) { print("using conservative Dmax")}
		lim.list<-vector('list', ncol(Combinations))
		for(j in 1:ncol(Combinations)){
			focal_grps<-groups[names(groups) %in% Combinations[,j]]
			focal_mrca<-lapply(focal_grps,function(x) {mmm<-getMRCA(phy,as.numeric(names(groups[groups == x])))
														if(is.null(mmm)) { if (VERBOSE == TRUE) {print("consult tree, some groups include only one taxon, using parent nodes to define limits...")}
														mmm <- max(nodepath(phy,as.numeric(names(groups[groups == x])),Ntip(phy)+ 1))}
														return(mmm)})
			focal_height<- min(edge[[1]][edge[[1]]$parent %in% unlist(focal_mrca),]$parent.height)
			if (VERBOSE == TRUE) {print("For focal tips:")
									print(phy$tip.label[as.numeric(names(focal_grps))])
									print("Dmax.t must occur before node height:")
									print(focal_height)
									print("______________________________________________________")}
			
			lim.list[[j]]<-focal_height
			}
		Combinations<-rbind(Combinations,unlist(lim.list))
		}
	}
	
	if (VERBOSE == TRUE)		{
		print("starting combinations...")
	}
	# calculate the C values 1-4 
	if(!is.null(groups) & conservative == TRUE){	## calculate Cmat under conservative Dmax option, passing oldest group mrca height to lim.height in calcCs
	Cmat <- apply(Combinations, 2, function(x) calcCsCt(tips = x[1:2], ancList=Ancs, allDists=allDists, phy=phy, VERBOSE=VERBOSE, allVals = allVals, edge = edge, lim.height = x[3]))
	Combinations<-Combinations[1:2,]
	} else{
	Cmat <- apply(Combinations, 2, function(x) calcCsCt(tips = x, ancList=Ancs, allDists=allDists, phy=phy, VERBOSE=VERBOSE, allVals = allVals, edge = edge))
	}	
	avgC <- apply(sapply(Cmat,"[[",2),1,mean)
	to_return<-list("mean" = avgC,"Cmat" = sapply(Cmat,"[[",2),"meas.path" = sapply(Cmat,"[",1))
	
	if(!is.null(groups)){
	grp_mat<-matrix(groups[match(Combinations,names(groups))],2,length(Combinations)/2)
	grp_combo<-combn(unique(groups),2)

	grp_result<-vector('list', ncol(grp_combo))
	grp_result[["overall"]]<-avgC
	for(m in 1:ncol(grp_combo)){
	focal_combo <- grp_combo[,m]
	focal_combo<-focal_combo[order(focal_combo)]
	focal_pairs<-which(apply(grp_mat,2,function(x) all(x[order(x)] == focal_combo)))
	grp_result[[paste(focal_combo,collapse = "_")]] <- apply(data.frame(sapply(Cmat,"[[",2)[,focal_pairs]),1,mean)
	}
	grp.temp<-do.call(cbind,grp_result)
	to_return[["grp.mean"]]<- cbind(grp.temp,overall.weighted = apply(grp.temp[,colnames(grp.temp) != "overall",drop = FALSE],1,mean))
	}
	return(to_return)
}
