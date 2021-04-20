\name{BASS.standard.generator}
\alias{BASS.standard.generator}
\title{Function that returns the results of standard Bass model}
\description{Function returning the results of a standard Bass model, BM, (Bass, 1969).
}

\usage{
BASS.standard.generator(m,p,q,tstart=1,n=50)
}

\arguments{
\item{m}{Market potential}
\item{p}{Coefficient of innovation}
\item{q}{Coefficient of imitation}
\item{tstart}{Inital time}
\item{n}{Number of observations}
}

\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}

\examples{
#Example 1

BASS.standard.generator(10,0.1,0.1,8,3)
}


\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.plot}}

\code{\link{make.instantaneous}}

\code{\link{BASS.generalized}}

\code{\link{GG.model}}

\code{\link{UCRCD}}

}

\keyword{BASS.standard.generator}

