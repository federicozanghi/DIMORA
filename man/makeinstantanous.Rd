\name{make.instantaneous}
\alias{make.instantaneous}
\title{Function that converts the data}
\description{Function that turns cumulative values into instantaneous values. This function may be useful in some cases, since all procedures in this package are built using instantaneous data.}

\usage{
make.instantaneous(cumulate.data)
}

\arguments{
\item{cumulate.data}{Cumulative data}
}

\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}

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

\code{\link{UCRCD}}

}

\keyword{BASS}

