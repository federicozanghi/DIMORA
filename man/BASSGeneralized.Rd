\name{BASS.generalized}
\alias{BASS.generalized}
\title{Function that estimates a Generalized Bass model}
\description{This function allows to estimate parameters of a Generalized Bass model, GBM, (Bass et al., 1994) with Exponential, Rectangular, Harmonic or Mixed shock, (Guseo et al. 2007).
}

\usage{
BASS.generalized(sales, shock = c("exp", "rett", "mixed", "harmonic"),
                 nshock, prelimestimates, alpha = 0.05, ous = 100,
                 display = T,max.iter=100)
}

\arguments{
\item{sales}{Instantaneous sales}
\item{shock}{Function characterizing the shock. Available options are: Exponential ('exp'), Rectangular ('rett'), Harmonic ('harmonic') and Mixed ('mixed').}
\item{nshock}{Number of shocks desired (from 1 to 3)}
\item{prelimestimates}{Vector containing the preliminary estimates of the parameters}
\item{alpha}{Desired significance level, the default value is 0.05}
\item{ous}{Numeric value for the out-of-sample forecasts, the default value is 100}
\item{display}{T or F, to display the plot of the model or not, the default value is T}
\item{max.iter}{Maximum number of iterations, the default value is 100}

}
\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}


\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{BASS.plot}}

\code{\link{make.instantaneous}}

\code{\link{GG.model}}

\code{\link{UCRCD}}

\code{\link{predict_SARMAX.refinement}}

}

\keyword{BASS.generalized}
