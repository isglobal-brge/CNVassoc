\name{getProbs}
\Rdversion{1.1}
\alias{getProbs}
\alias{getProbs.cghCall}
\alias{getProbs.cnv}


\title{
 Get posterior probabilities from an object of class 'cnv' or 'CGHcall'
}


\description{
 This function creates a list where each component correspond to a given probe. 
 Any component contains the posterior probabilities obtained from CGHcall algorithm.
}


\usage{
  getProbs(x)
  \method{getProbs}{cghCall}(x)
  \method{getProbs}{cnv}(x)
}

\arguments{
 \item{x}{an object of class 'CGHcall' or 'cnv'}
}


\value{
A list where each component correspond to a given probe. Any component contains the posterior probabilities obtained from CGHcall algorithm
}

\note{
See vignette for an example.
}

\author{
 This function was created using a script kindly provided by Mark van de Wiel
}


\keyword{utilities}

