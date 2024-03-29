% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onestart.R
\name{onestart}
\alias{onestart}
\title{Evaluate individual nodes for their value in invasion detection (invasion starting from a single node) in one simulation}
\usage{
onestart(adjmat, wtvec, start.choice, stoch = T)
}
\arguments{
\item{adjmat}{adjacency matrix for evaluation}

\item{wtvec}{vector of weights corresponding to nodes associated with the adjacency matrix adjmat}

\item{start.choice}{number of node where invasion starts}

\item{stoch}{logical var indicating whether adjacency matrix entries are fixed or probabilities}
}
\description{
For a given starting node (only one at this point), yields a vector of the number of nodes invaded for each sampling node by the time invasion reaches (is detected at) that sampling node, and the time step of first invasion at the sampling node - for one simulation.  Includes the option for differentially weighted probabilities of the existence of links.
}
\examples{
onestart(adjmat=Amat, wtvec= ___, start.choice=2) # assumes Amat exists - more complete example to be added
onestart(adjmat=sAmat, wtvec= ___, start.choice=2,stoch=T) # assumes Amat exists - more complete example to be added
}
\keyword{prioritization}
\keyword{sampling}
