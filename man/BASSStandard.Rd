\name{BASS.standard}
\alias{BASS.standard}
\title{Function that estimates a standard Bass model}
\description{This function allows to estimate the parameters of a standard Bass model, BM, (Bass, 1969) by using two different methods: nonlinear least squares, NLS, or OPTIM, which minimizes the residual sum of squares.
}


\usage{
BASS.standard(sales, method = "nls", prelimestimates = c(sum(sales) +
              100, 0.01, 0.1), ous= 100, alpha = 0.05, display = T, max.iter = 100)
}

\arguments{
\item{sales}{Instantaneous sales }
\item{method}{Method for parameter estimation, 'nls' or 'optim'}
\item{prelimestimates}{A vector containing the starting values for parameter estimation. Default values are m = sum(sales)+100, p = 0.01, q = 0.1}
\item{alpha}{Significance level for confidence intervals}
\item{display}{T or F, to display the plot of the model or not, the default value is T}
\item{ous}{Numeric value for the out-of-sample forecasts, the default value is 100}
\item{max.iter}{Maximum number of iterations, the default value is 100. This option works only with method="nls"}
}

\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}


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

\code{\link{UCRCD}}

}

\keyword{BASS.standard}

