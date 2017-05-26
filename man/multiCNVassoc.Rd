\name{multiCNVassoc}
\Rdversion{1.1}
\alias{multiCNVassoc}

\title{
Association between several CNVs and disease
}     

\description{
This function repeatedly calls CNVassoc function
}

\usage{
multiCNVassoc(x, formula, num.copies = 0:2, cnv.tol = 0.01, ...)
}

\arguments{
  \item{x}{a list of calling probabilities matrix for each CNV }
  \item{formula}{see 'formula' argument of \code{\link{CNVassoc}} function. }
  \item{num.copies}{See argument of \code{\link{cnv}} function. }
  \item{cnv.tol}{See argument of \code{\link{cnv}} function. }  
  \item{\dots}{other arguments passed through 'CNVassoc' function.}
}

\value{
A list of p-values for each CNV
}


\details{
See vignette for an example
}


\references{
Gonzalez JR, Subirana I, Escaramis G, Peraza S, Caceres A, Estivill X and Armengol L. Accounting for uncertainty when assessing association between copy number and disease: a latent class model. \emph{BMC Bioinformatics}, 2009;10:172.
}




\seealso{
\code{\link{CNVassoc}}
}


\keyword{misc}

