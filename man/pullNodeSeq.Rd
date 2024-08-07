% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pullNodeSeq.R
\name{pullNodeSeq}
\alias{pullNodeSeq}
\title{Extracts a vector of ancestors for a given taxon. Code written by Jonathan S. Mitchell for
Zelditch et al. (2017)}
\usage{
pullNodeSeq(phy, tip)
}
\arguments{
\item{phy}{The phylogeny of interest in phylo format}

\item{tip}{The tip of interest}
}
\value{
A vector of ancestors
}
\description{
pullNodeSeq Extracts a vector of all ancestors of a given taxon in a phylogeny.
}
\examples{

phylogeny<-rtree(100)
answer<-pullNodeSeq(phy=phylogeny,tip="t1")
}
\references{
Zelditch, M.L., J. Ye, J.S. Mitchell, and D.L. Swiderski. 2017. Rare ecomorphological
convergence on a complex adaptive landscape: Body size and diet mediate evolution of jaw shape in
squirrels (Sciuridae). Evolution 71:633-649.
}
