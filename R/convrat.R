#'Quantifies convergent evolution using the C1, C2, C3, and C4 measures as described by Stayton (2015).  
#'
#'Calculates the current phenotypic distance (Euclidean) between two taxa.  Then uses ancestral state reconstruction under a BM model to calculate the maximum phenotypic distance at any time between lineages leading from the most recent common ancestor of those two taxa to the tips.  Also calculate the total amount of phenotypic evolution in the clade defined by the most recent common ancestor of those lineages, and the total amount of phenotypic evolution in the input tree.  These quantities are used to calcualte C1-C4:  C1 = 1-(current distance / maximum ancestral distance); C2 = maximum ancestral distance - current distance; C3=C2/(total phenotypic evolution in the clade defined by the two taxa); C4 = C2/(total amount of phenotypic evolution in the entire tree).  If more than two convergent taxa are input, then C1-C4 are calculated for all possible pairs of taxa, and averaged.
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param convtips A list consisting of the names of all convergent taxa
#'
#'@details None
#'
#'@return Four numbers - C1, C2, C3, C4.   
#'
#'@import ape geiger MASS phytools 
#'
#'@export
#'
#'@references Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2013).
#'cluster: Cluster Analysis Basics and Extensions. R package version 1.14.4.
#'
#'Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics
#'and evolution in R langauge. Bioinformatics, 20, 289-290.
#'
#'Revell, L. J. (2012) phytools: An R package for phylogenetic comparative 
#'biology (and other things). Methods Ecol. Evol. 3 217-223.
#'
#'Stayton, C.T.  (2015).  The definition, recognition, and interpretation of
#'convergent evolution, and two new measure for quantifying and assessing the 
#'significance of convergence.  Evolution 69:2140-2153.
#'
#'@examples
#'
#'phyl<-rtree(10)
#'phendata<-fastBM(phyl,nsim=2)
#'convtips<-c("t1","t2","t3")
#'answer<-convrat(phyl,phendata,convtips)


convrat<-function(phyl,phendata,convtips)

{

#Error checking

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

if (nrow(phendata) != length(phyl$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    
#Different commands are required if there are only two putatively convergent 
#taxa, versus more than two.

ntax<-length(convtips)

if(ntax==2) {

	tipsdist<-sqrt(sum((phendata[convtips[1] ,]-phendata[convtips[2] ,])^2))

	mxdist<-maxdist(phyl,phendata,convtips[1],convtips[2])

	if (tipsdist>=mxdist) {
		C1<-0
		C2<-0
		}

	if (mxdist>tipsdist) {
		C1<-1-(tipsdist/mxdist)
		C2<-mxdist-tipsdist
		}

	

	lineagepaths<-ancestrallineages(phyl,phendata,convtips[1],convtips[2])

	t1totalevolution<-0
	t2totalevolution<-0

	for (i in 1:dim(lineagepaths[[1]])[1]-1) {

		branchdistance<-sqrt(sum(((lineagepaths[[1]][i,]-lineagepaths[[1]][i+1,])^2)))
		t1totalevolution<-t1totalevolution+branchdistance
		
		}

	for (i in 1:dim(lineagepaths[[2]])[1]-1) {

		branchdistance<-sqrt(sum(((lineagepaths[[2]][i,]-lineagepaths[[2]][i+1,])^2)))
		t2totalevolution<-t2totalevolution+branchdistance

		}

	totallineageevolution<-t1totalevolution+t2totalevolution

	C3<-C2/totallineageevolution
	
	cMRCA<-findMRCA(phyl,tips=convtips,type="node")
	subphyl<-extract.clade(phyl,cMRCA)	
	
	wholephylchanges<-sum(calcchanges(subphyl,phendata[subphyl$tip.label,]))

	C4<-C2/wholephylchanges

	}

if(ntax>2) {

	C1s<-c()
	C2s<-c()
	C3s<-c()
	C4s<-c()

	for (i in 1:ntax) {

		j<-i+1

		while (j<=ntax) {

			tipsdist<-sqrt(sum((phendata[convtips[i] ,]-phendata[convtips[j] ,])^2))

			mxdist<-maxdist(phyl,phendata,convtips[i],convtips[j])

			C1<-1-(tipsdist/mxdist)

			C2<-mxdist-tipsdist

			lineagepaths<-ancestrallineages(phyl,phendata,convtips[i],convtips[j])

			t1totalevolution<-0
			t2totalevolution<-0

			for (ii in 1:(dim(lineagepaths[[1]])[1]-1)) {

				branchdistance<-sqrt(sum(((lineagepaths[[1]][ii,]-lineagepaths[[1]][ii+1,])^2)))
				t1totalevolution<-t1totalevolution+branchdistance

				}

			for (ii in 1:dim(lineagepaths[[2]])[1]-1) {

				branchdistance<-sqrt(sum(((lineagepaths[[2]][ii,]-lineagepaths[[2]][ii+1,])^2)))
				t2totalevolution<-t2totalevolution+branchdistance

				}


			totallineageevolution<-t1totalevolution+t2totalevolution

			C3<-C2/totallineageevolution
	
			cMRCA<-findMRCA(phyl,tips=convtips,type="node")
			subphyl<-extract.clade(phyl,cMRCA)	
	
			wholephylchanges<-sum(calcchanges(subphyl,phendata[subphyl$tip.label,]))

			wholephylchanges<-sum(calcchanges(phyl,phendata))

			C4<-C2/wholephylchanges

			C1s<-c(C1s,C1)
			C2s<-c(C2s,C2)
			C3s<-c(C3s,C3)
			C4s<-c(C4s,C4)

			j<-j+1

			}

		}

	C1<-mean(C1s)
	C2<-mean(C2s)
	C3<-mean(C3s)
	C4<-mean(C4s)

	}

answer<-c(C1,C2,C3,C4)
names(answer)<-c("C1","C2","C3","C4")
answer

}