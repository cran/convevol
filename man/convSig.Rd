% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convSig.R
\name{convSig}
\alias{convSig}
\title{Uses simulations to assess the significance of C1-C4 measures of convergent
evolution as described in Stayton (2015). Code written by Jonathan S.
Mitchell for Zelditch et al. (2017)}
\usage{
convSig(phy, traits, focaltaxa, user.ace = NULL, nsim = 1000)
}
\arguments{
\item{phy}{The phylogeny of interest in phylo format}

\item{traits}{Phenotypic data for all tips}

\item{focaltaxa}{A list consisting of the names of all putatively convergent taxa}

\item{user.ace}{A matrix of user supplied ancestral trait values at internal nodes (formatted as "traits" but with node number as rownames)}

\item{nsim}{The number of simulations to use to assess significance}
}
\value{
C1-C4 convergence measures for all pairs of putatively convergent
taxa and their associated p-values.
}
\description{
convSig  calculates the significance of measures of convergent evolution.
}
\details{
This script simulates data according to a Brownian motion model of
evolution, and then assesses convergene on that simulated data. The number
of times that the simulated data produces greater convergence than that
observed in the empirical data is used to calculate a P-value.
}
\examples{

phy<-rtree(100)
traits<-fastBM(phy,nsim=3)
focaltaxa<-c("t1","t50","t100")
answer<-convSig(phy,traits,focaltaxa,user.ace=NULL,nsim=10)
}
\references{
Stayton, C.T. 2015. The definition, recognition, and interpretation of convergent evolution, and
two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140-2453.

Zelditch, M.L., J. Ye, J.S. Mitchell, and D.L. Swiderski. 2017. Rare ecomorphological
convergence on a complex adaptive landscape: Body size and diet mediate evolution of jaw shape in
squirrels (Sciuridae). Evolution 71:633-649.
}
