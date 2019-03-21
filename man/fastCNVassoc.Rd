\name{fastCNVassoc}
\alias{fastCNVassoc}
                              
\title{Fast association analysis between a CNV or imputed SNPs and phenotype for case-control and cohort GWA studies.}

\description{
  This function performs an association analyses with several genetic variants with uncertainty (CNVs or imputed SNPs) and a response, maybe adjusting for covariates (e.g. clinical covariates, stratification, ...). It uses the Newthon-Raphson procedure with analytic likelihood derivatives and it has written in C language to speed up the process making feasible to analyse hundreds of thousands of variants, It also incorporates the possibility to use several cores to make calculations even faster.
}


\usage{
  fastCNVassoc(probs, formula, data, model = "additive",
    family = "binomial", nclass = 3, colskip = 5, tol = 1e-06, max.iter = 30, 
    verbose = FALSE, multicores=0)
}

\arguments{
  \item{probs}{either a matrix containing the probabilities of genetic variants in IMPUTE format (i.e. each rows represent a variant and every 'nclass' columns an individual. The first 'colskip' columns refers to the variant info (e.g. position, rs name, alleles, etc.).}  
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. In the right side of ~ covariates must be included. If no covariates are present in the model just type 1.}
  \item{data}{an optional data frame, list or environment (or object coercible by 'as.data.frame' to a data frame) containing the variables in the model.  If not found in 'data', the variables are taken from 'environment(formula)'.}
  \item{model}{Genetic model to be tested. Possible values are "multiplicative" (model free, e.g. co-dominant) or "additive", partial matching allowed. Only "additive" model is implemented.}
  \item{family}{a description of the error distribution and link function to be used in the model. This must be a character string naming a family function. Possible values are "binomial" and "weibull". Default value is "binomial"}
  \item{nclass}{integer specifying the number of possibles alleles or genotypes. Default value is 3 (tipycally for SNPs).}
  \item{colskip }{integer specifying the number of columns to be skipped in order to read the probabilities. This columns may contain the SNPs info (such as rs name, position, alleles and chromosome). Default value is 5 for IMPUTE format files.}
  \item{tol}{Tolerance for convergence in fitting model. Default value is 1e-06.}
  \item{max.iter}{Maximum number of iterations in fitting model. Default value is 30.}
  \item{verbose}{logical. If TRUE the number of current analysed variants is shown in the console. Default value is FALSE}
  \item{multicores}{integer indicating the number of cores to be used. It uses 'parallel'. Default value is 0 indicating that only one core is used and 'parallel' package is not required. For Windows OS, 'multicores'>1 is not supported.}
}



\value{  

 A data.frame with the following variables:
 
  - variant: consecutive integer from one to the number of analyzed variants (CNV or imputed SNPs) in the same order as in the probs matrix.

  - beta coefficient: log-Odds Ratio for binary response or log-Hazard Ratio for time-to-event response.

  - se: standard error of beta coefficient.

  - zscore: ratio between beta and se

  - pvalue: p-value of association between each genetic variant and response, maybe adjusted by covariates.

  - iter: number of iterations necessary achieve convergence during the Newton-Raphson algorithm used to fit the model.
  
  See examples for further illustration about all previous issues. 

}


\note{

  The order of individuals from probabilities ('probs' argument) matrix must be the same as in the response and covariates variables.
  
  The 'subset' and 'na.action' is not implemented. Therefore, no missings are allowed in probabilities, response or covariates.
  
  It is important be aware whether the number of iterations has achieved the maximum (30 by default). In this case, the results may be not reliable.
  
  The probability matrix ('probs') must have 'nclass' * N + 'colskip' columns, where N is the number of individuals.

}


\references{

Subirana I, Gonz?lez JR. Genetic Association Analysis and Meta-Analysis of Imputed SNPs in Longitudinal Studies. \emph{Genet Epidemiol}, 2013 Jul;37(5):465-77.   \cr

Gonzalez JR, Subirana I, Escaramis G, Peraza S, Caceres A, Estivill X and Armengol L. Accounting for uncertainty when assessing association between copy number and disease: a latent class model. \emph{BMC Bioinformatics}, 2009;10:172. 

}



\seealso{

\code{\link{CNVassoc}}, \code{\link{multiCNVassoc}}, \code{\link{CNVtest}}

}


\examples{

\dontrun{

require(CNVassoc)
require(parallel)

# read imputed SNP probabilities from a file-.
# Example from SNPTEST software of 500 cases and 500 controls on 200 imputed SNPS.
fileprobs <- system.file("exdata/SNPTEST.probs",package="CNVassoc")

# build response (500 controls and 500 cases).
resp<-rep(0:1,each=500)

# generate two covariates randombly
N<-1000
covar1<-rnorm(N)  # contiuous covariate
covar2<-factor(sample(1:3,N,replace=TRUE),labels=c("A","B","C")) # categorical covariate

# run with 6 cores. Under Windows OS, multicore must be <=1.
system.time(
res<-fastCNVassoc(fileprobs,resp~covar1+covar2,family="binomial",multicore=6)
)
res

# build a time-to-event response randomly
set.seed(123456)
times <- rexp(N,1)
cens <- rbinom(N,1,0.8)

system.time(
res<-fastCNVassoc(fileprobs,Surv(times,cens)~covar1+covar2,family="weibull",multicore=6)
)
res


}


}

\keyword{misc}

