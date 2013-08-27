\name{convnum}
\alias{convnum}
\title{Quantify convergence by the number of convergent events}
\usage{
  convnum(phyl, phendata, convtips)
}
\arguments{
  \item{phyl}{The phylogeny of interest in phylo format}

  \item{phendata}{Phenotypic data for all tips}

  \item{convtips}{A list consisting of the names of all
  convergent taxa}
}
\value{
  The number of lineages that have crossed into the region
  of trait space occupied by the convergent taxa.
}
\description{
  This program takes in a set of taxa that are already
  suspected to be convergent in a particular area of
  morphospace.  It then counts the number of times that a
  lineage has invaded that region of morphospace.
  Essentially soubles the sum of two numbers
}
\details{
  This function will construct an ellipse around all
  convergent taxa.  Then it will reconstruct ancestral
  states throughout the phylogeny, and use those to
  determine how many lineages have crossed into this
  ellipse from the outside.
}
\examples{
phyl<-rtree(10)
phendata<-fastBM(phyl,nsim=2)
convtips<-c("t1","t2","t3")
answer<-convnum(phyl,phendata,convtips)
}
\references{
  Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M.,
  Hornik, K.(2013). cluster: Cluster Analysis Basics and
  Extensions. R package version 1.14.4.

  Paradis, E., J. Claude, and K. Strimmer (2004) APE:
  Analyses of phylogenetics and evolution in R langauge.
  Bioinformatics, 20, 289-290.

  Revell, L. J. (2012) phytools: An R package for
  phylogenetic comparative biology (and other things).
  Methods Ecol. Evol. 3 217-223.
}

