\name{simCNVdataWeibull}
\alias{simCNVdataWeibull}

\title{
Simulate of CNV and a right censored Weibull distributed trait
}
\description{
This function simulates intensity for a CNV and a time to event response (followed-up cohort study design) for different scenarios 
}
\usage{
simCNVdataWeibull(n, mu.surrog, sd.surrog, w, lambda, shape, time.cens = Inf, 
                  cnv.random = FALSE)
}

\arguments{
  \item{n}{An integer indicating the desired number of individuals to be simulated}
  \item{mu.surrog}{A vector containing the signal (surrogate variable) means for every copy number status (latent classes). Its length must be equal to the number of latent classes}
  \item{sd.surrog}{A vector containing the signal standard deviation for every copy number status. Its length must be equal to mu.surrog.}
  \item{w}{A vector containing the frequencies for every copy number status. Its length must be equal to mu.surrog and its components must sum up one.}
  \item{lambda}{A vector containing the means of the response variable for every copy number status}
  \item{shape}{A vector containing the shape of the response variable for every copy number status}
  \item{time.cens}{Censoring time, e.g. end of follow-up}
  \item{cnv.random}{A logical value. TRUE means that copy number status is drawn under a multinomial distribution with proportions indicated by 'w'. FALSE means that the real simulated frequency is always the same and is rounded to the most similar integer to the frequencies indicated by 'w'. Default value is FALSE}
}

\value{
  Data frame with individual simulated data per row and with the following variables:
  \item{resp}{Time to event or censoring variable (response) }
  \item{cens}{Censoring indicator} 
  \item{surrog}{Signal intensity following a mixture of normals with means, standard deviations
     and proportions specified by \code{mu.surrog}, \code{sd.surrog} and \code{w} respectively}
  \item{cnv}{True copy number status}
}


\seealso{
  \code{\link{simCNVdataBinary}}, 
  \code{\link{simCNVdataCaseCon}},
  \code{\link{simCNVdataPois}},
  \code{\link{simCNVdataNorm}},
  \code{\link{cnv}},
  \code{\link{CNVassoc}}
}
\examples{
library(survival)  
maf<-0.3
hr<-1.5
set.seed(123)
simData<-simCNVdataWeibull(n=3000, mu.surrog=c(0,0.5,1), sd.surrog=rep(0.15,3), 
         w=c((1-maf)^2,2*maf*(1-maf), maf^2), lambda=0.05*c(1,hr,hr^2), shape=rep(1,3), 
         time.cens=1.5, cnv.random = FALSE)
CNV<-cnv(simData$surrog,mix.method="EMmixt")
getQualityScore(CNV,type="CNVtools")  
mod<-CNVassoc(Surv(resp, cens)~CNV,data=simData,family="weibull")
CNVtest(mod)
summary(mod)
}



\keyword{utilities}
