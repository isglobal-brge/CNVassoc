\name{dataMLPA}
\docType{data}
\alias{dataMLPA}
\title{MLPA data}

\usage{data(dataMLPA)}

\description{
  This data set contains data from a MLPA assay for a case-control study. The data has intensities for two genes that can be used to infer copy number status. It also contains the TRUE copy number status obtained by using PCR technology. 
}


\format{
  dataMLPA is a data.frame with the following columns:
  \tabular{ll}{
	id \tab The unique identifiers of individuals \cr
	casco \tab Case-control stauts 0:control 1:case \cr
	Gene1 \tab Intensities for Gene1 \cr
	Gene2 \tab Intensities for Gene2 \cr
	PCR.Gene1 \tab True copy number status for Gene1\cr
	PCR.Gene2 \tab True copy number status for Gene2\cr
	quanti \tab Simulated continuos variable. It was generated to illustrate how to perform association analysis between a CNV and continuous traits \cr
	cov \tab Simulated continuous variable. It was generated to illustrate how to perform association adjusting for covariates \cr
  }

}


\references{
Gonzalez JR, Subirana I, Escaramis G, Peraza S, Caceres A, Estivill X and Armengol L. Accounting for uncertainty when assessing association between copy number and disease: a latent class model. \emph{BMC Bioinformatics}, 2009;10:172.
}

\source{
  This data set was generated and kindly provided by Estivill's lab. The data is still unpublished and it has been made available in a blinded format for reproducing the findings presented in the paper 
}

\keyword{datasets}

