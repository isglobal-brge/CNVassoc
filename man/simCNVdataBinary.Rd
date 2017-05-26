\name{simCNVdataBinary}
\Rdversion{1.1}
\alias{simCNVdataBinary}

\title{
Simulation of CNV and discrete traits
}

\description{
This function simulates intensity for a CNV and a binary trait response for different scenarios
}

\usage{
simCNVdataBinary(n, mu.surrog, sd.surrog, w, p0, or, cnv.random = FALSE)
}

\arguments{
  \item{n}{number of simulated individuals}
  \item{mu.surrog}{a vector of intensity signal means for every copy number status}
  \item{sd.surrog}{a vector of intensity signal standard deviations for every copy number status}
  \item{w}{a vector of copy number status proportions}
  \item{p0}{prevalence of disease (trait) for populations with zero copies (reference category)}
  \item{or}{a vector of odds ratio for one, two,... copies respect to zero copies}
  \item{cnv.random}{A logical value. TRUE means that copy number status is drawn under a multinomial distribution with proportions indicated by 'w'. FALSE means that the real simulated frequency is always the same and is rounded to the most similar integer to the frequencies indicated by 'w'. Default value is FALSE}
}


\value{
  Data frame with individual simulated data per row and with the following variables:
  \item{resp}{Trait (response) variable following a Bernoulli distribution given the CNV status}
  \item{surrog}{Signal intensity following a mixture of normals with means, standard deviations
     and proportions specified by \code{mu.surrog}, \code{sd.surrog} and \code{w} respectively.}
  \item{cnv}{True copy number status}
}

\details{
  This function is useful to calculate the power of association models with binary traits under different scenarios, e.g. setting different degrees of association (odds ratios), considering different degrees of uncertainty controlled by the distribution of intensity signal data, i.e. mean \code{mu.surrog}, standard deviation \code{sd.surrog} and proportion \code{w}, etc.
}


\seealso{
  \code{\link{simCNVdataCaseCon}}, 
  \code{\link{simCNVdataNorm}},
  \code{\link{simCNVdataPois}},
  \code{\link{simCNVdataWeibull}},
  \code{\link{cnv}},
  \code{\link{CNVassoc}}
}

\examples{
maf<-0.3
set.seed(123)
simData<-simCNVdataBinary(n=1000, mu.surrog=c(0,0.5,1), sd.surrog=rep(0.15,3), 
         w=c((1-maf)^2,2*maf*(1-maf),maf^2), p0=0.1, or=c(1.3,1.3^2), cnv.random = FALSE)
CNV<-cnv(simData$surrog,mix.method="EMmixt")
getQualityScore(CNV,type="CNVtools")
mod<-CNVassoc(resp~CNV,data=simData,family="binomial")
CNVtest(mod)
summary(mod)
}

\keyword{utilities}

                         