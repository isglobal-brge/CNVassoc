\name{NeveData}
\docType{data}
\alias{NeveData}
\title{Breast Cancer aCGH data}
\usage{data(NeveData)}
\description{
This data set contains breast cancer data studied by Neve et al. (2006). The data consists on CGH arrays of 1MB resolution. The authors chose the 50 samples that could be matched to the name tokens of caArrayDB data (June 9th 2007).

This object is a list with two components. The first component corresponds to a data.frame containing 2621 rows and 54 columns with aCGH data (4 columns for the annotation and 50 log2ratio intensities). The second component contains information about strogen receptor positivity (dichotomous variable; 0: negative, 1: positive)
}


\format{
  The first component is a data.frame with the following columns:
  \tabular{ll}{
	Clone \tab The unique identifiers of array elements \cr
	Chrom \tab Chromosome number of each array element \cr
	kb \tab Chromosomal position in bp of each array element \cr
	kb.1 \tab Chromosomal position in bp of each array element \cr
	X600MPE \tab Raw log2 ratios for breast cancer sample X600MPE \cr
	AU565 \tab Raw log2 ratios for breast cancer sample AU565 \cr
	  ...        \tab  ...      \cr
        ZR75B \tab Raw log2 ratios for breast cancer sample ZR75B \cr
  }

  The second component is a vector with information about strogen receptor positivity
}

\source{
 Data are freely available from the bioconductor website (http://www.bioconductor.org/) through the package Neve2006
}

\keyword{datasets}

