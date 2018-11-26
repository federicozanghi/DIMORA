\name{BASS.standard}
\alias{BASS.standard}
\title{Function that returns estimate parameters of bass standard model and the plot}
\description{This function allows you to calculate the parameters of the bass model in 2 different ways:
using nls or using optim(which minimize the residual squared sum).
}


\usage{
BASS.standard(sales,method,prelimestimates,alpha,display)
}

\arguments{
\item{sales}{The istantanuous sales}
\item{method}{The way you decide to estimate parameters, 'nlm' or 'optim'}
\item{prelimestimates}{A vector containing the starting values, that the algorithm use to estimate the parameters, the default values are: m = sum(sales)+100, p = 0.01, q = 0.1 }
\item{alpha}{The significativity level for confidence intervals}
\item{display}{T or F if you want to display the graphs of the models or not, the default value is T}
}

\author{
Zanghi Federico \email{federico.zanghi@studenti.unipd.it}
Department of Statistics, University of Padova, Italy \url{http://www.unipd.it}
}



\examples{
#Example 1
Sales <- c(840,1470,2110,4000,7590,10950,10530,9470,7790,5890)
BASS.standard(sales = Sales,method = 'nls')
}

\seealso{

\code{\link{BASS.standard.generator}}

\code{\link{BASS.plot}}

\code{\link{make.instantaneous}}

\code{\link{BASS.generalized}}

\code{\link{GG.model}}

}

\keyword{BASS.standard}

