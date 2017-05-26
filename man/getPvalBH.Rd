\name{getPvalBH}
\Rdversion{1.1}
\alias{getPvalBH}

\title{
 Corrected p values using Benjamini & Hochberg approach
}
\description{
 This functions corrects the association p-values using the Benjamini & Hochberg approachby for multiple testing.
}

\usage{
getPvalBH(x)
}


\arguments{
  \item{x}{a list containing p values}
}

\details{
This function calls 'p.adjust' to compute 'BH' correction
}

\value{
A data frame with the blocks and corrected p-values
}



\seealso{
\code{\link[stats]{p.adjust}}
}


\keyword{misc}

