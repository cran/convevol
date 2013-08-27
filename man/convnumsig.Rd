\name{convnumsig}
\alias{convnumsig}
\title{Assess the significance of convergent evolution using simulations and the convnum metric}
\usage{
  convnumsig(phyl, phendata, convtips, nsim)
}
\arguments{
  \item{phyl}{The phylogeny of interest in phylo format}

  \item{phendata}{Phenotypic data for all tips}

  \item{convtips}{A list consisting of the names of all
  convergent taxa}

  \item{nsim}{The number of simulatons to conduct}
}
\value{
  A list, consisting first of the p-value for the observed
  convnum, and second of a vector containing all of the
  simulated convnum values.  Also displays a histogram of
  all of the simulated convnum values.
}
\description{
  Simulates evolution along a given phylogeny, using
  parameters derived from observed data, and calculates the
  convnum metric for each simulation for a set of
  user-defined taxa.  Then compares the observed convnum
  value to the simulated values to assess the significance
  of the observed levels of convergent evolution.
}
\details{
  None
}
\examples{
phyl<-rtree(10)
phendata<-fastBM(phyl,nsim=2)
convtips<-c("t1","t2","t3")
answer<-convnumsig(phyl,phendata,convtips,10)
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

