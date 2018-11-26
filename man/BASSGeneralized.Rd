\name{BASS.generalized}
\alias{BASS.generalized}
\title{Function that returns the results of generalized Bass model}
\description{This function allows you to estimate the generalized Bass model with an Esponential, Rectangular, Armonic or Mixed shock.
}

\usage{
BASS.generalized(sales,shock,nshock,prelimestimates,alpha)
}

\arguments{
\item{sales}{The istantanuous sales}
\item{shock}{The function which caracterize the shock, the implemented options are: Esponential ('exp'), Rectangular ('rett'), Armonic ('armonic') and Mixed ('mixed').}
\item{nshock}{The number of shocks desired (from 1 to 3) }
\item{prelimestimates}{A vector containing the preliminar estimates of the parameters}
\item{alpha}{The desired significativity, the defaul value is 0.05}

}
\author{
Zanghi Federico \email{federico.zanghi@studenti.unipd.it}
Department of Statistics, University of Padova, Italy \url{http://www.unipd.it}
}


\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{BASS.plot}}

\code{\link{make.instantaneous}}

\code{\link{GG.model}}

}

\keyword{BASS.generalized}
