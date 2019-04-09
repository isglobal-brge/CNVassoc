\name{SNPTEST}
\alias{SNPTEST}
\alias{cases}
\alias{controls}
\title{Case-control data with SNPTEST format}
\usage{data(SNPTEST)}
\description{
Case-control data for 200 individuals (100 cases and 100 controls) and 50 SNPs 
using SNPTEST format (e.g. probabilities of AA, AB and BB genotypes)
}
\format{
Two objects: cases and controls. Each object is a matrix with 100 rows and 1505 
columns. Column 1-5 correspond to annotation. Files 6-1505 correspond to the 
probabilities of being AA, AB or BB for 500 SNPs
}
\source{
\url{https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html}
}
\value{
Two data.frames with cases and controls. Each data.frame consists of 50 rows 
(individuals, and 1505 columns with imputed probabilities, and SNP info in the 
first five columns).
}
\examples{
data(SNPTEST)
cases[1:10, 1:10]
controls[1:10, 1:10]
}
\keyword{datasets}
