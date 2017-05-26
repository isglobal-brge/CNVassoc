\name{plotSignal}
\alias{plotSignal}
\title{
plots the intensities of a CNV univariate signal data
}
\description{
This function creates a plot with probe intensity
}
\usage{
plotSignal(x, my.colors = c("black", "red", "blue"), ylab = "Peak Intensity", 
      xlab = c("individuals", "Phenotype"), case.control = NULL, cex.leg = 0.6, 
      dens.bw = "nrd0", dens.adjust = 1, n = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector with probe intensities}
  \item{my.colors}{Colours for each copy number status.}
  \item{ylab}{Label of y-axis}
  \item{xlab}{Label of x-axis}
  \item{case.control}{Vector indicating case-control status}
  \item{cex.leg}{Size of legend}
  \item{dens.bw}{Adjustment for intensity signal density curve. See argument 'bw' of \code{\link[stats]{density}} function for more details}
  \item{dens.adjust}{Adjustment for intensity signal density curve. See argument 'adjust' of \code{\link[stats]{density}} function for more details}
  \item{n}{integer indicating the number of points to be placed on the plot interactively (using \code{locator})) to define the thresholds that separate the different copy number status, colouring the points differently according to the assigned copy number status. If it zero or negative it makes no possible to place the threshold points.}
  \item{\dots}{Other arguments passed to \code{plot.default}}
}

\details{
See vignette for further description
}

\references{
Gonzalez JR, Subirana I, Escaramis G, Peraza S, Caceres A, Estivill X and Armengol L. Accounting for uncertainty when assessing association between copy number and disease: a latent class model. \emph{BMC Bioinformatics}, 2009;10:172.
}

\seealso{
\code{\link{cnv}}
}


\examples{
data(dataMLPA)
plotSignal(dataMLPA$Gene2)
}


\keyword{utilities}
