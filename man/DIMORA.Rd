\name{DIMORA}
\alias{DIMORA}
\title{DIMORA Package}
\description{
In the following package some innovation diffusion models have been implemented to analyze nonlinear growth processes, for both univariate and bivariate time series data.

The implemented models are: standard Bass model, Generalized Bass model (with rectangular shock, exponential shock, mixed shock and harmonic shock, 1 to 3 shocks available), Dynamic market potential model, and UCRCD model.

The Bass model, BM, (Bass, 1969), consists of a simple differential equation that describes the process of how new products get adopted in a population, while the Generalized Bass model, GBM, (Bass et al., 1994) is a generalization of the Bass model with a function x(t), capturing the changing speed of diffusion.
In some real processes the market potential may be not constant over time and a dynamic market potential model is needed. The Guseo-Guidolin model, GGM, (Guseo and Guidolin, 2009) is a specification of this situation.

The UCRCD model, (Guseo and Mortarino, 2014) is a diffusion model useful to capture the dynamics of competition between two products within the same market.
See the full documentation for more details.}





\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}

Andrea Savio \email{svandr97@gmail.com}
}

\examples{
library(DIMORA)
}

\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{BASS.plot}}

\code{\link{BASS.generalized}}

\code{\link{GG.model}}

\code{\link{UCRCD}}

\code{\link{make.instantaneous}}

\code{\link{SARMAX.refinement}}

\code{\link{classic.plot.sarima}}

\code{\link{summary.Dimora}}

\code{\link{predict.Dimora}}

\code{\link{plot.Dimora}}



}

\keyword{DIMORA}

