\name{CNVassoc}
\Rdversion{1.1}
\alias{CNVassoc}
\alias{anova.CNVassoc}
\alias{logLik.CNVassoc}
\alias{print.CNVassoc}
\alias{summary.CNVassoc}
\alias{print.anova.CNVassoc}
\alias{print.summary.CNVassoc}
                               
\title{Association analysis between a CNV and phenotype}

\description{
  This function performs an association analysis between a CNV and a dependent variable (phenotype) using a latent class model that incorporates the uncertainty arising from calling procedure. The phenotype may be quantitative or categorical. In the second case (e.g. case-control studies) this variable must be coded as 1 (for cases) and 0 (for controls). The association can be adjusted for other covariates (e.g. clinical covariates, stratification, ...) 
}


\usage{
  CNVassoc(formula, data, subset, na.action, model = "multiplicative", 
           family = "binomial", tol = 1e-06, max.iter = 30, emsteps = 0, 
           verbose = FALSE, coef.start, sigma.start, alpha.start=1)
}

\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Right side of ~ should have an object of class 'cnv'.}
  \item{data}{an optional data frame, list or environment (or object coercible by 'as.data.frame' to a data frame) containing the variables in the model.  If not found in 'data', the variables are taken from 'environment(formula)'.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain 'NA's.  The default is set by the 'na.action'        setting of 'options', and is 'na.fail' if that is unset.  The 'factory-fresh' default is 'na.omit'. Another possible value is 'NULL', no action.  Value
          'na.exclude' can be useful.} 
  \item{model}{Genetic model to be tested. Possible values are "multiplicative" (model free, e.g. co-dominant) or "additive", partial matching allowed. Default value is "multiplicative".}
  \item{family}{a description of the error distribution and link function to
          be used in the model. This must be a character string naming a
          family function. Possible values are "binomial", "gaussian", "poisson" or "weibull". Default value is "binomial"}
  \item{tol}{Tolerance for convergence in fitting model. Default value is 1e-06.}
  \item{max.iter}{Maximum number of iterations in fitting model. Default value is 30.}
  \item{emsteps}{Number of iterations using Expectation Maximization (EM) alghorithm to set initial values before using Newton-Rapson (NR) in fitting model. Default value is zero, that means that EM step is not performed}
  \item{verbose}{logical. If TRUE parameter values for each iteration are shown in the console. Default value is FALSE}
  \item{coef.start}{initial values for coefficients in NR procedure}
  \item{sigma.start}{initial values for scale parameter (only for "gaussian") in NR procedure}
  \item{alpha.start}{initial values for shape parameter (only for "weibull") in NR procedure}
}



\value{
 An object of class 'CNVassoc'. 

  'print' returns model parameter estimates

  'summary' returns a summary table similar to summary.glm 

  'anova' performs a Likelihood Ratio Test comparing two nested models fitted using \code{CNVassoc}

  'logLik' returns the log-likelihood of a model fitted using \code{CNVassoc}
  
  See examples for further illustration about all previous issues. 
}



\references{
Gonzalez JR, Subirana I, Escaramis G, Peraza S, Caceres A, Estivill X and Armengol L. Accounting for uncertainty when assessing association between copy number and disease: a latent class model. \emph{BMC Bioinformatics}, 2009;10:172. 
}




\seealso{
\code{\link{cnv}},\code{\link{CNVtest}}
}

\examples{
data(dataMLPA)
CNV  <-  cnv(x  =  dataMLPA$Gene2,  threshold.0  =  0.01,  mix.method  =  "mixdist")
modmul  <-  CNVassoc(casco  ~  CNV,  data  =  dataMLPA,  model  =  "mul")
modmul
summary(modmul)
anova(modmul,update(modmul,model="add"))
logLik(modmul)
}

\keyword{misc}

