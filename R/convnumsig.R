#'Assess the significance of convergent evolution using simulations and the convnum metric
#'
#'Simulates evolution along a given phylogeny, using parameters derived from observed data, and calculates the convnum metric for each simulation for a set of user-defined taxa.  Then compares the observed convnum value to the simulated values to assess the significance of the observed levels of convergent evolution.  
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param convtips A list consisting of the names of all convergent taxa
#'@param nsim The number of simulatons to conduct
#'
#'@details None
#'
#'@return A list, consisting first of the p-value for the observed convnum, and second of a vector containing all of the simulated convnum values.  Also displays a histogram of all of the simulated convnum values.
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
#'@examples
#'
#'phyl<-rtree(10)
#'phendata<-fastBM(phyl,nsim=2)
#'convtips<-c("t1","t2","t3")
#'answer<-convnumsig(phyl,phendata,convtips,10)

convnumsig<-function(phyl,phendata,convtips,nsim)

#Evaluates the significance of convergence as measured by convnum.  

{

#Then the observed value, and the ellipse.

phendata<-as.matrix(phendata)

convtaxa<-phendata[unlist(convtips) ,]

convell<-ellipsoidhull(convtaxa,0.0001)

ob<-convnum(phyl,phendata,convtips)

#Then the simulations

ancvals<-multianc(phyl,phendata)

rootvals<-ancvals[length(phyl$tip.label) ,]

C<-vcv.phylo(phyl)

vcv<-phyl.vcv(phendata,C,0)

simdata<-sim.char(phyl,vcv$R,nsim,model=c("BM"),root=rootvals)

#And then all of the assessments of those simulations

sobs<-c()
moresig<-0

nbran<-dim(phyl$edge)

for (i in 1:nsim) {

	sphendata<-simdata[, , i]

	salldata<-multianc(phyl,sphendata)

	#phylomorphospace(phyl,sphendata)

	#plotellipse(convell)

	cross<-0

	j<-1

	while (j<=nbran[1]) {

		isinanc=FALSE
		isindes=FALSE
	
		anc<-phyl$edge[j,1]
		des<-phyl$edge[j,2]

		ancval<-salldata[anc ,]
		desval<-salldata[des ,]

		ancdev<-ancval-convell$loc
		desdev<-desval-convell$loc

		cutoff<-1.00*convell$d2  #cutoff<-(1+tol)*convell$d2

		if (t(ancdev)%*%ginv(convell$cov)%*%ancdev<=cutoff) {isinanc=TRUE}else{isinanc=FALSE}

		if (t(desdev)%*%ginv(convell$cov)%*%desdev<=cutoff) {isindes=TRUE}else{isindes=FALSE}

		if (isinanc!=TRUE & isindes==TRUE) {

			cross<-cross+1

			#And we'll highlight those branches in red.

			#pts<-rbind(ancval[1:2],desval[1:2])

			#lines(pts,col="red")

			#arrows(ancval[1],ancval[2],desval[1],desval[2],length=0.1,angle=30,code=2,col="red")

			}

		j<-j+1

		}
	
	if(cross>=ob) {moresig<-moresig+1}
	
	sobs<-c(sobs,cross)

	}

hist(sobs)

p<-moresig/nsim

all<-list(p,sobs)
}
		