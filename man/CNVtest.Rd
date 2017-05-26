\name{CNVtest}
\Rdversion{1.1}
\alias{CNVtest}
\alias{print.CNVtest}

\title{
 Testing association between a CNV and phenotype
}

\description{
  This function perform a Wald test or a Likelihood Ratio Test (LRT) to determine whether a CNV is associated with the phenotype.
}

\usage{
  CNVtest(x, type = "Wald")
  \method{print}{CNVtest}(x, ...)
}

\arguments{
  \item{x}{An object of class 'CNVassoc'}
  \item{type}{The statistical test used. Possible values are "Wald" for Wald test or "LRT" for Likelihood Ratio Test}
  \item{\ldots}{Further arguments passed to or from other methods}
}

\value{
  An object of class 'CNVtest' with methods for the generic function \code{print}, returning the Chi-squared value, its degrees of freedom and the corresponding p-value of CNV significance in the associated test fitted by \code{CNVassoc} function.
}

\references{
Gonzalez JR, Subirana I, Escaramis G, Peraza S, Caceres A, Estivill X and Armengol L. Accounting for uncertainty when assessing association between copy number and disease: a latent class model. \emph{BMC Bioinformatics}, 2009;10:172. 
}


\seealso{
\code{\link{CNVassoc}},\code{\link{cnv}}
}

\examples{
data(dataMLPA)
CNV<- cnv(x = dataMLPA$Gene2, threshold.0 = 0.01, mix.method = "mixdist")
mod<-CNVassoc(formula = casco ~ CNV, data = dataMLPA, model = "mul")
CNVtest(mod, type = "LRT")
CNVtest(mod, type = "Wald")
}

\keyword{misc}

