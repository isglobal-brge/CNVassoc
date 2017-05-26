\name{simCNVdataNorm}
\Rdversion{1.1}
\alias{simCNVdataNorm}


\title{
Simulation of CNV and quantitative traits
}

\description{
This function simulates intensity for a CNV and a quantitative trait response for different scenarios
}

\usage{
simCNVdataNorm(n, mu.surrog, sd.surrog, w, mu.y, sd.y, cnv.random = FALSE)
}

\arguments{
  \item{n}{An integer indicating the desired number of individuals to be simulated}
  \item{mu.surrog}{A vector containing the signal (surrogate variable) means for every copy number status (latent classes). Its length must be equal to the number of latent classes}
  \item{sd.surrog}{A vector containing the signal standard deviation for every copy number status. Its length must be equal to mu.surrog.}
  \item{w}{A vector containing the frequencies for every copy number status. Its length must be equal to mu.surrog and its components must sum up one.}
  \item{mu.y}{A vector containing the means of the response variable for every copy number status. Its length must be equal to mu.surrog.}
  \item{sd.y}{A single number indicating the residual standard deviation}
  \item{cnv.random}{A logical value. TRUE means that copy number status is drawn under a multinomial distribution with proportions indicated by 'w'. FALSE means that the real simulated frequency is always the same and is rounded to the most similar integer to the frequencies indicated by 'w'. Default value is FALSE}
}

\value{
  Data frame with individual simulated data per row and with the following variables:
  \item{resp}{Continous trait variable (response) }
  \item{surrog}{Signal intensity following a mixture of normals with means, standard deviations
     and proportions specified by \code{mu.surrog}, \code{sd.surrog} and \code{w} respectively}
  \item{cnv}{True copy number status}
}

\details{
  This function is useful to calculate the power of association models for a continuous (normal-distributed) trait under different scenarios ,e.g. setting different degrees of association (effects), considering different degrees of uncertainty controlled by the distribution of intensity signal data, i.e. mean \code{mu.surrog}, standard deviation \code{sd.surrog} and proportion \code{w}, etc.
}


\seealso{
  \code{\link{simCNVdataBinary}}, 
  \code{\link{simCNVdataCaseCon}},
  \code{\link{simCNVdataPois}},
  \code{\link{simCNVdataWeibull}},
  \code{\link{cnv}},
  \code{\link{CNVassoc}}
}
\examples{
set.seed(123)
maf<-0.3
effect<-3
simData<-simCNVdataNorm(n=1000, mu.surrog=c(0,0.5,1), sd.surrog=rep(0.15,3), 
         w=c((1-maf)^2,2*maf*(1-maf), maf^2), mu.y=100+c(0,effect,2*effect), 
         sd.y=rep(20,3), cnv.random = FALSE)
CNV<-cnv(simData$surrog,mix.method="EMmixt")
getQualityScore(CNV,type="CNVtools")
mod<-CNVassoc(resp~CNV,data=simData,family="gaussian",emsteps=10)
CNVtest(mod)
summary(mod)
}

\keyword{misc}

