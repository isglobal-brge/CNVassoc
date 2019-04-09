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
See vignette for more examples example
}
\references{
Gonzalez JR, Subirana I, Escaramis G, Peraza S, Caceres A, Estivill X and 
Armengol L. Accounting for uncertainty when assessing association between 
copy number and disease: a latent class model. \emph{BMC Bioinformatics}, 
2009;10:172.
}
\seealso{
\code{\link{CNVassoc}}
}
\examples{
data(SNPTEST)
nSNP <- nrow(cases)
probs <- lapply(1:nSNP, function(i) {
    snpi.cases <- matrix(as.double(cases[i, 6:ncol(cases)]), ncol = 3, 
                        byrow = TRUE)
    snpi.controls <- matrix(as.double(controls[i, 6:ncol(controls)]), ncol = 3, 
                        byrow = TRUE)
    return(rbind(snpi.cases, snpi.controls))
})
casecon <- rep(1:0, c(500, 500))
pvals <- multiCNVassoc(probs, formula = "casecon~CNV", model = "add",
                        num.copies = 0:2, cnv.tol = 0.001)
pvalsBH <- getPvalBH(pvals)
head(pvalsBH)
}
\keyword{misc}
