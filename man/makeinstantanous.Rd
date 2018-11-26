\name{make.instantaneous}
\alias{make.instantaneous}
\title{Function that transforms your data}
\description{Function that transforms your cumulate data into instantaneous, this function can be usefull because all the function in this package reqire the istantanus data.}

\usage{
make.instantaneous(cumulate.data)
}

\arguments{
\item{cumulate.data}{ The cumulate data}
}

\author{
Zanghi Federico \email{federico.zanghi@studenti.unipd.it}
Department of Statistics, University of Padova, Italy \url{http://www.unipd.it}
}


\examples{
#Example 1
cumulate.data = c(1,2,3,6,12)
data.inst = make.instantaneous(cumulate.data)
}

\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{BASS.plot}}

\code{\link{BASS.generalized}}

\code{\link{GG.model}}

}

\keyword{BASS}

