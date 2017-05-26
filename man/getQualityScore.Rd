\name{getQualityScore}
\Rdversion{1.1}
\alias{getQualityScore}
\alias{getQualityScore.default}
\alias{getQualityScore.cnv}
\alias{print.QualityScore}

\title{
Computes a quality score for a CNV fit
}

\description{
 This function provides different types of  measurements of uncertainty after CNV calling
}

\usage{
getQualityScore(x, ...)
\method{getQualityScore}{default}(x, sds, w, type, iter = 10000, threshold = 0.1, ...)
\method{getQualityScore}{cnv}(x, type = "class", iter = 10000, threshold = 0.1, ...)
}

\arguments{
  \item{x}{and object of class \code{cnv} or means vector of intensity signal for each copy number status}
  \item{...}{further arguments passed to or from \code{getQualityScore} methods}
  \item{type}{the type of quality score measurement computed. Possible values are "class", "CNVtools" or "CANARY" (see Details)}
  \item{iter}{number of iterations when \code{type=''class''} or \code{type=''CANARY''} is specified}
  \item{threshold}{a value to compute the proportion of sample individuals with confidence score bigger 
     than it (see Details)}
  \item{sds}{standard deviations vector of intensity signal for each copy number status}
  \item{w}{copy number status proportions vector}
}


\details{
 The quality scores measures how well the clusters are separated. 
 If compares the locations of the means with the standard error for each pair of adjacent cluster. 
 Obviously, except for probability of good classification (\code{type="class")}, 
 the lower quality score the highest uncertainty.
 There are 3 possible types of quality score measurements:
 "class": probability of good classification),
 "CNVtools": the score defined in 'CNVtools' package) and 
 "CANARY": proportions of sample individuals with confidence score bigger than \code{threshold}. The confidence
 score is defined as the ratio between the second biggest copy number call probability divided by
 the biggest one.
}


\value{
An object of class \code{getQualityScore} with a single number of quality score.
}




\note{
For \code{cnv} objects created directly from probabilities and not from fitting a univariate intensity signal,
only "class" quality score type can be calculated.
}



\examples{
data(dataMLPA)
CNV<-cnv(x = dataMLPA$Gene2, threshold.0 = 0.01, mix.method = "mixdist")
getQualityScore(CNV,type="class")
getQualityScore(CNV,type="CNVtools")
getQualityScore(CNV,type="CANARY")
}

\keyword{ misc}

