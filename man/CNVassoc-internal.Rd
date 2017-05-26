% --- CNVassoc-internal.Rd ---

\name{CNVassoc-internal}
\alias{CNVassoc-internal}
\alias{assocCGH.i}
\alias{EMlogistic}
\alias{EMnorm}
\alias{hessianLinear}
\alias{hessianLogistic}
\alias{linear.fit}
\alias{logistic.fit}
\alias{mix}
\alias{NRlogistic}
\alias{NRnorm}
\alias{plot.cnv.intensities}
\alias{plot.cnv.probabilities}
\alias{matrix2vector}
\alias{vector2matrix}
\alias{is.quantitative}
\alias{bicMix}
\alias{plapply}
\alias{mixture}

\title{Internal CNVassoc functions}
\description{Internal CNVassoc functions}
\usage{
assocCNV.i(x, formula, num.copies, cnv.tol, ...)
EMlogistic(y, X, w, beta, tol = 10^-6, max.iter = 1000, verbose = FALSE)
EMnorm(y, X, w, beta, sigma, variant, tol = 10^-6, max.iter = 1000, verbose = FALSE)
EMpoisson(y, X, w, beta, tol = 10^-6, max.iter = 1000, verbose = FALSE)
EMWeibull(y, cens, X, w, beta, alpha, tol = 10^-6, max.iter = 1000, verbose = FALSE)
hessianLinear(beta, sigma, y, w, X, variant)
hessianLogistic(beta, y, w, X, variant)
hessianPoisson(beta, y, w, X, variant)
hessianWeibull(beta, alpha, y, cens, w, X, variant)
ifelsem(test, yes, no)
linear.fit(x, y, weights, tol = 1e-08, max.iter = 25, verbose = FALSE)
logistic.fit(x, y, weights, tol = 1e-08, max.iter = 25, verbose = FALSE)
matrix2vector(betam, variant)
mix(mixdat, method, num.class, mean.ini, sigma.ini, pi.ini, var.equal)
mixture(intensities, num.class, mix.method, threshold.0, threshold.k, mu.ini,
        sigma.ini, pi.ini, var.equal)
NRlogistic(y, X, w, beta, variant, tol, max.iter, verbose)
NRnorm(y, X, w, beta, sigma, variant, tol = 10^-6, max.iter = 1000, verbose = FALSE)
plot.cnv.intensities(x, my.colors = c("black", "red", "blue"), ylab = "Peak Intensity", 
        xlab = c("individuals", "Phenotype"), case.control, cex.leg = 0.8, 
        dens.bw = "nrd0", dens.adjust = 1, ...)
plot.cnv.probabilities(x, my.colors = c("black", "red", "blue"), case.control, 
         ylab = "CNV probability", xlab, ...)
vector2matrix(betav, variant, J)
is.quantitative(formula, data)
getProbsRegions.i(i, blocks, probs, annotation, nclass)
plapply(X, FUN, ...)
translate(obj)
}

\details{These are not to be called by the user}
\keyword{internal}
