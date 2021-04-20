\name{classic.plot.sarima}
\alias{classic.plot.sarima}
\title{Plot for \code{SARMAX.refinement} function.}
\description{Function to obtain the same graphs that are obtained with "display = T" in the SARMAX.refinement function.
}

\usage{
classic.plot.sarima(model,typeog=0,...)
}

\arguments{
\item{model}{An object of class \code{Dimora}, i.e. \code{BASS.generalized}, \code{BASS.standard} or \code{GG.model}.}
\item{typeog}{Three plots available, residuals versus time (\code{typeog = 1)}, residual ACF (\code{typeog = 2}), fitted values on real data (\code{typeog = 3}), mode \code{typeog = 0} shows all plots. The default value is \code{typeog = 0}.}
\item{...}{Other graphic parameters.}
}
\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}

\keyword{classic.plot.sarima}
