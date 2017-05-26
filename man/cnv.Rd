\name{cnv}
\Rdversion{1.1}
\alias{cnv}
\alias{cnvBatches}
\alias{cnvDefault}
\alias{is.cnv}
\alias{plot.cnv}
\alias{print.cnv}

\title{CNV object}

\description{
  \code{cnv} creates a 'cnv'  object
  \code{is} returns \code{TRUE} if \code{x} is of class 'cnv'
  \code{print} gives a summary for an object of class 'cnv' including ...
  \code{plot}  plots an object of class 'cnv' ...
 
}

\usage{
  cnv(x, batches, ...)
  cnvDefault(x, num.copies, num.class, cnv.tol = 0.001, mix.method = "mixdist", 
             check.probs = TRUE,  threshold.0, threshold.k, mu.ini, sigma.ini, pi.ini, 
             cutoffs = NULL, check.alpha = 0.05, check.cnv = TRUE, var.equal)
  cnvBatches(intensities, batches, threshold.0, threshold.k, common.pi = TRUE, ...)
  is.cnv(obj)
  \method{plot}{cnv}(x, ...)
  \method{print}{cnv}(x, digits = 4, ...)
}


\arguments{
  \item{x}{a vector of CNV intensity signal for each individual, or a matrix with CNV calling probabilities per row}
  \item{num.copies}{vector with copy number status values, i.e, number of copies or a vector of characters indicating loss ('l'), normal ('n') or gain ('g') for example} 
  \item{num.class}{integer indicating how many classes CNV contains} 
  \item{cnv.tol}{error tolerance when x is a probability matrix and row sums are not identical to one}
  \item{mix.method}{normal mixture fitting method when x is a vector of univariate CNV signal intensities. Current methods are "mixdist" that uses the function \code{mix} from the package \code{mixdist}, "mclust" that uses de function \code{Mclust} from the package \code{mclust} and "EMmixt" that uses an internal function \code{EMmixt} from the \code{CNVassoc} package. The last two are based on Expectation-Maximization procedure and the first one is based on quasi-Newton-Raphson procedure}
  \item{check.probs}{logical. If TRUE it checks weather row sums are equal to one +/- \code{cnv.tol} when x is a probability matrix}  
  \item{threshold.0}{assigns zero copies (or first copy number status) to all individuals whose CNV signal intensity is lower than threshold.0} 
  \item{threshold.k}{assigns k copies (or last copy number status) to all individuals whose CNV signal intensity is bigger than threshold.k} 
  \item{mu.ini}{an opcional vector to specify the initial values of means when fitting a normal mixture to CNV intensity signal data}
  \item{sigma.ini}{an opcional vector to specify the initial values of standard deviations when fitting a normal mixture to CNV intensity signal data} 
  \item{pi.ini}{an optional vector to specify the initial values of copy number status probabilities when fitting a normal mixture to CNV intensity signal data} 
  \item{cutoffs}{a vector indicating the cut-off points to assign the copy number status assign individuals to the individuals according to the categories defined by these cut-off points on CNV intensity signal data}
  \item{check.alpha}{significance level to goodness-of-fit test indicating weather the normal mixture model to CNV intensity data has been fitted appropriately}
  \item{check.cnv}{logical. If TRUE, cnv functions returns and error when normal mixture model does not fit well to the univariate CNV intensity signal data} 
  \item{var.equal}{logical. If TRUE, standard deviation are supposed to be the same for all copy number status when fitting univariate CNV intensity signal data}
  \item{intensities}{a vector with the univariate CNV intensity signal data}
  \item{batches}{a vector indicating the batch (leave it missing if no batch effect is present)}
  \item{common.pi}{logical. If TRUE, copy number status probabilities for each individual are computed estimating specific means and standard deviations separately for every batch, but the same population copy number status probabilities for all batches. It is suggested to leave it as TRUE}
  \item{obj}{an object of any class}
  \item{digits}{number of digits when printing a \code{cnv} object}
  \item{\ldots}{other arguments passed to \code{cnvDefault}, \code{print.default} or \code{plot.cnv}. The arguments passed to \code{plot.cnv}are the same as the ones for the \code{\link{plotSignal}} function}
}


\details{
When argument batches is not specified, then cnvDefault is used, otherwise cnvBatch is called.
If univariate CNV intensity signal data is used to create the cnv class object, then one can introduce the batch effect if it necessary. But, if other algorithms have been used previously and the cnv class object is created directly from the CNV calling probabilities matrix, then it is not possible to specify the batch  argument.
The batch effect is important when cases and controls have been genotyped in different platforms for example. In this situations, the platform should be introduced in the batch argument as a vector indicating which platform every CNV intensity signal data comes from.
Generic \code{plot} function applied on a 'cnv' class object performs two types of plots whether 'cnv' class object has been created from univariate CNV intensity signal data or whether it has been created directly from a probability matrix provided by any CNV calling algorithm. The first type is a plot similar to the one created by \code{\link{plotSignal}} function, and the second type is a barplot.
}


\value{
\code{cnv} return an object of class 'cnv' with generic function such as \code{print} or \code{plot} implemented for this kind of objects.
\code{is.cnv} is a function that returns TRUE of FALSE weather \code{obj} is of class 'cnv' or not.
}


\references{
Gonzalez JR, Subirana I, Escaramis G, Peraza S, Caceres A, Estivill X and Armengol L. Accounting for uncertainty when assessing association between copy number and disease: a latent class model. \emph{BMC Bioinformatics}, 2009;10:172. 
}


\seealso{
\code{\link{CNVassoc}}, \code{\link{plotSignal}}
}

\examples{
data(dataMLPA)
CNV  <-  cnv(x  =  dataMLPA$Gene2,  threshold.0  =  0.01,  mix.method  =  "mixdist")
CNV
plot(CNV)
}

\keyword{misc}

