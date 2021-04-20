\name{SARMAX.refinement}
\alias{SARMAX.refinement}
\title{SARMAX refinement of a \code{Dimora} model.}
\description{This function allows to refine results obtained with a model of the \code{Dimora} class through a SARIMAX model, (Guseo et. al., 2007). The function \code{predict_SARMAX.refinement} is aviable for this procedure.}

\usage{
SARMAX.refinement(fit,arima_order,seasonal_order,period=1,display=F)
}

\arguments{
\item{fit}{Object of class\code{Dimora}, i.e.   \code{BASS.generalized}, code{BASS.standard} or code{GG.model}.}
\item{arima_order}{Non-seasonal part of the SARIMA model }
\item{seasonal_order}{Seasonal part of the SARIMA model}
\item{period}{Period of the seasonal part of the ARIMA model. Default value is 1 (daily data).}
\item{display}{Logical "T" and "F", "T" for some graphs to be shown}
}

\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}


\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{BASS.plot}}

\code{\link{make.instantaneous}}

\code{\link{GG.model}}

\code{\link{predict_SARMAX.refinement}}

\code{\link{UCRCD}}
}

\keyword{SARMAX.refinement}
