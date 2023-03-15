#'Computes Ct values for a pair of tips. Internal, called in calcConv.
#'
#'calcCsCt Computes Ct values for a pair of tips. Internal, called in calcConv.
#'
#'@param tips vector of two tips
#'@param ancList list of node paths for all tips in the user provided phylogeny
#'@param allDists matrix of phenotypic distances between all exterior and interior node pairs
#'@param phy The time calibrated phylogeny of interest in phylo format
#'@param VERBOSE logical value indicating whether model information should be printed during computation
#'@param allVals a matrix of observed and reconstructed phenotypes for all user supplied traits at interior and exterior nodes
#'@param edge a list of data frames, each including the edge matrix of user supplied phylogeny, along with node heights and reconstructed phenotype for each of the user supplied traits 
#'@param lim.height an optional tree height used to limit Dmax.t, passed only if groups are defined and a conservative test is run (see calcConv)
#'
#'@details Function incorporates the optimizations introduced by Zelditch et al. (2017), which significantly improve runtimes
#'@details Reconstructions part way along branches are obtained using equation [2] of Felsenstein (1985), following code modified from the phytools (Revell, 2012) function contMap
#'
#'@return A list of the Ct values
#'
#'@import MASS phytools
#'
#'@importFrom stats na.omit
#'
#'@export
#'
#'@references Grossnickle DM, Brightly WH, Weaver LN, Stanchak KE, Roston RA, Pevsner SK, Stayton CT, Polly PD, Law CJ. 2022. A cautionary note on quantitative measures of phenotypic convergence. in revision
#'Zelditch ML, Ye J, Mitchell JS, Swiderski DL. 2017. Rare ecomorphological convergence on a complex adaptive landscape: Body size and diet mediate evolution of jaw shape in squirrels (Sciuridae). Evolution 71: 633-649
#'Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution and two new measures for quantifying and assessing the significance of convergence. Evolution 69(8): 2140-2153.
#'Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.
#'Felsenstein, J. 1985. Phylogenies and the comparative method. American Naturalist, 125, 1-15.
#'

calcCsCt <- function(tips, ancList, allDists, phy, VERBOSE=FALSE, allVals, edge, lim.height = NULL)	{
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
	# delete - ancPairs <- sapply(unit1, function(t1) sapply(unit2, function(t2) allDists[t1,t2]))
	
	# save edge info and create meas.path object with node and interpolated branch reconstructed values
	edge_parent <- rownames(edge[[1]])[edge[[1]][,1] %in% c(unit1,unit2)]
	edge_child  <- rownames(edge[[1]])[edge[[1]][,2] %in% c(unit1,unit2)]
	edge_path<-as.numeric(intersect(edge_parent,edge_child))
	path_info<-lapply(edge,function(x) x[edge_path,])
	meas.path<-rbind(data.frame(node = unit1,path="P1"),
					data.frame(node = unit2,path="P2"))
	meas.path[meas.path$node %in% c(tips[1],tips[2]),]$path <-"tip"
	meas.path[meas.path$node == MRCA,]$path <- "mrca"
	meas.path$height <- sapply(1:nrow(meas.path),function(x) nodeheight(phy,meas.path$node[x]))
	meas.path<- cbind(meas.path,allVals[meas.path$node,])
	colnames(meas.path)[4:ncol(meas.path)]<-gsub("^","node.anc",1:ncol(allVals))	## revise to preserve original trait names?
	meas.path<-unique(meas.path)
	interp<-vector('list', nrow(meas.path))
	for(i in 1:nrow(meas.path)){
		if(!meas.path$path[i] %in% c("tip","mrca")){
		focal.node<-meas.path[i,]$node
		focal.height<-meas.path[i,]$height
	
		temp<-lapply(path_info,function(x) x[!x[["parent"]] == focal.node & !x[["child"]] == focal.node & x[["parent.height"]] < focal.height & x[["child.height"]] > focal.height,])
		if(unique(unlist(lapply(temp,nrow))) %in% 0){interp[[i]] <- rep(NA,ncol(allVals))
		}else{	
		
		interp[[i]]<- unlist(lapply(temp,function(x) {
				A1<-x["parent.anc"]
				A2<-x["child.anc"]
				t1<-focal.height-x["parent.height"]
				b<-x["child.height"] - x["parent.height"]
	
				return((A1/t1 + A2/(b-t1))/(1/t1 + 1/(b-t1)))
				}))
		
		}
		} else {interp[[i]] <- rep(NA,ncol(allVals))}
	}
	meas.path<-cbind(meas.path,do.call(rbind, interp))
	colnames(meas.path)[(4+ncol(allVals)):ncol(meas.path)]<-gsub("^","int.anc",1:ncol(allVals))
	
	diff_list<-vector('list', nrow(meas.path))
	for(i in 1:nrow(meas.path)){
	temp<-meas.path[i,]
		if(!temp$path %in% c("tip","mrca")){
		diff_list[[i]]<-dist(rbind(unlist(temp[,4:(4+ncol(allVals)-1)]),unlist(temp[,(4+ncol(allVals)):ncol(meas.path)])))
		} else(diff_list[[i]]<-NA)
	}
	meas.path$anc.diff<-unlist(diff_list)
	meas.path[meas.path$path == "tip",]$anc.diff<-Dtip
	
	# sum changes along each lineage including the MRCA
	lin1 <- c(unit1, as.character(MRCA))
	Lin1 <- sum(sapply(1:(length(lin1)-1), function(x) allDists[lin1[x],lin1[x+1]]))
	lin2 <- c(unit2, as.character(MRCA))
	Lin2 <- sum(sapply(1:(length(lin2)-1), function(x) allDists[lin2[x],lin2[x+1]]))
	
	# sum all pairwise distances between nodes of the subtree from the MRCA of tip1 & tip2 (deprecated, see new version below)
	##totalClade <- as.character(getDescendants(phy, MRCA))
	##totalMove <- sum(allDists[totalClade, totalClade]) / 2
	
	## calculate new version of totalMove based off of evolutionary distance covered over each branch of subtree from the MRCA of tip1 & tip2
	focal.edge<-which(edge[[1]]$child %in% getDescendants(phy,MRCA))
	C4.edge<-lapply(edge,function(x) x[focal.edge,])
	rows<-nrow(C4.edge[[1]])
	ptraits<-sapply(1:rows,function(x) lapply(edge,function(y) y[x,"parent.anc"]))
	ctraits<-sapply(1:rows,function(x) lapply(edge,function(y) y[x,"child.anc"]))
	totalMove<-sum(sapply(1:rows,function(z) dist(rbind(ptraits[,z],ctraits[,z]))))
	
	# Dmax as the maximum distance between ancestor pairs
	
	##Below includes two options for treating Dmax, either allowing it to equal Dtip or preventing this. If preventing, results in Inf for C2 - C4
	##whenever sister taxa are included, note too the ability to measure divergence is also limited by the fact that Dmax in diverging taxa will (almost?)
	##always be between the most recent internal nodes, and thus measure will merely reflect branch length (likely result of sampling intensity).
	
	##Dmax <- unique(max(abs(na.omit(meas.path$anc.diff))))	## allow Dmax = Dtip
	
	if(!is.null(lim.height)) {
	Dmax <- max(na.omit(meas.path[!meas.path$path == "tip" & meas.path$height < lim.height,]$anc.diff))	## option for limiting Dmax age if conservative measure is being used
	} else {
	Dmax <- max(na.omit(meas.path[!meas.path$path == "tip",]$anc.diff))	## disallow Dmax = Dtip
	}
		
	C1 <- 1 - (Dtip / Dmax)
	C2 <- Dmax - Dtip
	C3 <- C2 / (Lin1 + Lin2)
	C4 <- C2 / totalMove
	return(list(meas.path,c(Ct1=C1, Ct2=C2, Ct3=C3, Ct4=C4)))
}
