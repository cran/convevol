#'Extracts lineages leading to two tips, t1 and t2, from their most recent common ancestor.
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param t1 The first tip of interest
#'@param t2 The second tip of interest
#'
#'@details None
#'
#'@return A list containing two matrices.  Each matrix corresponds to a tip.  The matrix 
#'consists of reconstructed ancestral values for all nodes leading from the mrca of both
#'tips to the tip.      
#'
#'@import ape geiger phytools 
#'
#'@export
#'
#'@references Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics
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
#'answer<-ancestrallineages(phyl,phendata,"t1","t2")


ancestrallineages<-function(phyl,phendata,t1,t2)


{

#Error checking

if (!inherits(phyl, "phylo")) {
        stop("your tree must be class 'phylo'")
    }

if (nrow(phendata) != length(phyl$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    
if (is.null(rownames(phendata))) {
	warning("no row names for data.  Assuming that the rows are in the same order as tips.")
      rownames(X) <- phyl$tip.label
	}

if (is.finite(t1)) {}
	else {t1<-labelstonumbers(phyl,t1)}

if (is.finite(t2)) {}
	else {t2<-labelstonumbers(phyl,t2)}

if (t1>length(phyl$tip))
	stop("your first tip isn't in the phylogeny.")

if (t2>length(phyl$tip))
	stop("your second tip isn't in the phylogeny.")




#The function

alldata<-multianc(phyl,phendata)

anctimes<-node.depth.edgelength(phyl)

combineddata<-cbind(anctimes,alldata)

#Then we go through and grab just the mrca of our two tips, and all nodes 
#between the mrca and the tips.  This part will depend on findanc and will
#be very similar to what is done in maxdist.

mrcas<-mrca(phyl)

mrcat1t2<-mrcas[t1,t2]

#First, the part for t1

t1path<-combineddata[t1 ,]

anc<-findanc(phyl,t1)

t1path<-rbind(t1path, combineddata[anc[1] ,])

while (anc[1] != mrcat1t2) {

	anc<-findanc(phyl, anc[1])

	t1path<-rbind(t1path, combineddata[anc[1] ,])

	}

#And then something identical for t2

t2path<-combineddata[t2 ,]

anc<-findanc(phyl,t2)

t2path<-rbind(t2path, combineddata[anc[1] ,])

while (anc[1] != mrcat1t2) {

	anc<-findanc(phyl, anc[1])

	t2path<-rbind(t2path, combineddata[anc[1] ,])

	}

#We just want the trait values for now

dims<-dim(t1path)

t1pathtraits<-t1path[, 2:dims[2]]
t2pathtraits<-t2path[, 2:dims[2]]

#And we'll bind them together

answer<-list(t1pathtraits,t2pathtraits)

}
