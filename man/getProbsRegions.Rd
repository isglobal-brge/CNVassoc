\name{getProbsRegions}
\Rdversion{1.1}
\alias{getProbsRegions}


\title{
Get posterior probabilities for blocks/regions 
}

\description{
See vignette for further details.
}

\usage{
getProbsRegions(probs,regions,intensities,nclass=3)
}

\arguments{
  \item{probs}{probabilities from CGH calling algorithms}
  \item{regions}{probs that define a segment}
  \item{intensities}{mean probe intensities for each region}
  \item{nclass}{See 'nclass' argument of function \code{\link[CGHcall]{CGHcall}}}
}


\note{
See vignette for an example.
}


\keyword{misc}

