\name{BASS.standard.generator}
\alias{BASS.standard.generator}
\title{Function that returns the results of standard Bass model}
\description{Function that returns the results of standard Bass model.
}

\usage{
BASS.standard.generator(m,p,q,tstart=1,n=50)
}

\arguments{
\item{m}{potential of market}
\item{p}{parameter of inovations}
\item{q}{parameter of imitators}
\item{tstart}{inital time}
\item{n}{number of values}
}

\author{
Zanghi Federico \email{federico.zanghi@studenti.unipd.it}
Department of Statistics, University of Padova, Italy \url{http://www.unipd.it}
}

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

}

\keyword{BASS.standard.generator}

