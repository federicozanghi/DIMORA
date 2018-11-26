\name{BASS.plot}
\alias{BASS.plot}
\title{Function that draws the plot of Bass model}
\description{Function that draws the plot of Bass model.
}

\usage{
BASS.plot(data)
}

\arguments{
\item{data}{result of BASS.standard.generator}
}

\author{
Zanghi Federico \email{federico.zanghi@studenti.unipd.it}
Department of Statistics, University of Padova, Italy \url{http://www.unipd.it}
}


\examples{
#Example 1

data = BASS.standard.generator(10,0.1,0.1,8,3)
BASS.plot(data)
}

\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{make.instantaneous}}

\code{\link{BASS.generalized}}

\code{\link{GG.model}}

}

\keyword{BASS}

