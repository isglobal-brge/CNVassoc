\name{getPvalBH}
\Rdversion{1.1}
\alias{getPvalBH}
\title{
    Corrected p values using Benjamini & Hochberg approach
}
\description{
    This functions corrects the association p-values using the Benjamini & 
    Hochberg approachby for multiple testing.
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
