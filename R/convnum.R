#'Quantify convergence by the number of convergent events
#'
#'This program takes in a set of taxa that are already suspected to be convergent in a particular area of  morphospace.  It then counts the number of times that a lineage has invaded that region of morphospace.  Essentially soubles the sum of two numbers
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param convtips A list consisting of the names of all convergent taxa
#'
#'@details This function will construct an ellipse around all convergent taxa.  Then it will reconstruct ancestral states throughout the phylogeny, and use those to determine how many lineages have crossed into this ellipse from the outside.  
#'
#'@return The number of lineages that have crossed into the region of trait space occupied by the convergent taxa.  
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
#'answer<-convnum(phyl,phendata,convtips)

convnum <-
function(phyl,phendata,convtips)

#convnum is a different approach to convergence.  Rather than looking at 
#relative distances between tips or lineages, this program takes in a set of
#taxa that are already suspected to be convergent in a particular area of 
#morphospace.  It then counts the number of times that a lineage has invaded
#that region of morphospace.

{

if (length(unlist(convtips))<=ncol(phendata)) {
stop("You must have fewer variables than putatively convergent taxa")
}

#First, defining the area of morphospace by an ellipse.

phendata<-as.matrix(phendata)

convtaxa<-phendata[unlist(convtips) ,]

convell<-ellipsoidhull(convtaxa,0.0001)
plotell<-ellipsoidhull(convtaxa[, 1:2],0.0001)

#And then getting the ancestral states

#dims<-dim(phendata)

#allancstates<-matrix(data=0, nr=dims[1]-1, nc=dims[2])

i<-1

while (i<=ncol(phendata)) {

tempdata<-phendata[, i]

tempancstates<-fastAnc(phyl,tempdata)

if(i==1){allancstates<-matrix(0,length(tempancstates),ncol(phendata))}

allancstates[, i]<-tempancstates

i<-i+1

}

alldata<-rbind(phendata,allancstates)

#Draw a phylomorphospace, which'll come in handy later.

phyl$node.label<-NULL

phylomorphospace(phyl,phendata[, 1:2],label=FALSE)

points(phendata[, 1:2],col="yellow")

points(phendata[unlist(convtips), 1:2],col="red")

#plotellipse(plotell)

#Now to go through all of the lineages

cross<-0

nbran<-dim(phyl$edge)

crossedges<-c()

i<-1

while (i<=nbran[1]) {

isinanc=FALSE
isindes=FALSE

anc<-phyl$edge[i,1]
des<-phyl$edge[i,2]

ancval<-alldata[anc ,]
desval<-alldata[des ,]

ancdev<-ancval-convell$loc
desdev<-desval-convell$loc

cutoff<-1.00*convell$d2  #cutoff<-(1+tol)*convell$d2

if (t(ancdev)%*%ginv(convell$cov)%*%ancdev<=cutoff) {isinanc=TRUE}else{isinanc=FALSE}

if (t(desdev)%*%ginv(convell$cov)%*%desdev<=cutoff) {isindes=TRUE}else{isindes=FALSE}

if (isinanc!=TRUE & isindes==TRUE) {

cross<-cross+1

#crossedges<-c(crossedges,i)

#And we'll highlight those branches in red.

#pts<-rbind(ancval[1:2],desval[1:2])

#lines(pts,col="red")

arrows(ancval[1],ancval[2],desval[1],desval[2],length=0.1,angle=30,code=2,col="red")

}

i<-i+1

}

cross

}
