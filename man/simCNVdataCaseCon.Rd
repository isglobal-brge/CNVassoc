\name{simCNVdataCaseCon}
\alias{simCNVdataCaseCon}

\title{
Simulation of CNV in a case-control study design
}
\description{
This function simulates intensity for a CNV within cases and control groups for different scenarios
}
\usage{
simCNVdataCaseCon(n0, n1, w0, or, mu.surrog0, sd.surrog0, mu.surrog1 = mu.surrog0, 
    sd.surrog1 = sd.surrog0, random = TRUE)
}


\arguments{
  \item{n0}{number of controls simulated}
  \item{n1}{number of cases simulated}
  \item{w0}{vector of proportions of copy number status in controls}
  \item{or}{a vector of odds ratio for one, two,... copies respect to zero copies}
  \item{mu.surrog0}{vector of means of CNV intensity signal, per copy number status, in control group}
  \item{sd.surrog0}{vector of standard deviations of CNV intensity signal, per copy number status, in control group}
  \item{mu.surrog1}{vector of means of CNV intensity signal, per copy number status, in control group}
  \item{sd.surrog1}{vector of standard deviations of CNV intensity signal, per copy number status, in control group}
  \item{random}{A logical value. TRUE means that individuals (rows) are randomly permuted, and FALSE means that simulated 'data.frame' contains controls first and then cases. Default value is TRUE}
}

\value{
  Data frame with individual simulated data per row and with the following variables:
  \item{resp}{Trait (response) variable with 0 or 1 if the individual is a control or a case
      respectively}
  \item{surrog}{Signal intensity following a mixture of normals with means, standard deviations
     and proportions specified by \code{mu.surrog}, \code{sd.surrog} and \code{w} respectively, 
     within cases and controls}
  \item{cnv}{True copy number status}
}

\details{
  This function is useful to calculate the power of association models in a case control study design under different scenarios             ,e.g. setting different degrees of association (odds ratios), considering different degrees of uncertainty controlled by the distribution of intensity signal data, i.e. mean \code{mu.surrog}, standard deviation \code{sd.surrog} and proportion \code{w}, etc.
}

\seealso{
  \code{\link{simCNVdataBinary}}, 
  \code{\link{simCNVdataNorm}},
  \code{\link{simCNVdataPois}},
  \code{\link{simCNVdataWeibull}},
  \code{\link{cnv}},
  \code{\link{CNVassoc}}
}
\examples{
maf<-0.3
set.seed(123)
simData<-simCNVdataCaseCon(n0=1000, n1=1000, mu.surrog0=c(0,0.5,1), sd.surrog0=rep(0.15,3), 
         mu.surrog1=c(0,0.5,1), sd.surrog1=rep(0.15,3), 
         w0=c((1-maf)^2,2*maf*(1-maf), maf^2), or=c(1.3,1.3^2),  
         random = FALSE)
CNV<-cnv(simData$surrog,mix.method="EMmixt")
getQualityScore(CNV,type="CNVtools")
mod<-CNVassoc(resp~CNV,data=simData,family="binomial")
CNVtest(mod)
summary(mod)
}

\keyword{utilities}
