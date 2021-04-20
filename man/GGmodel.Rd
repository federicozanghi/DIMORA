\name{GG.model}
\alias{GG.model}
\title{Function that estimates a dynamic market potential model }
\description{This function allows to estimate the parameters of a model with a time-dependent market potential, m(t).  The market potential may be defined according to the form proposed in Guseo and Guidolin (2009), generating the GG model, GGM. Other forms for m(t) may be defined according to the following: m(t) must depend on t and be a cdf.

To use the function, two options are available:

- Use the classic GGM, by setting function with sales, preliminary estimates and alpha

- Use the model with another m(t) function, by setting function with sales, preliminary estimates, mt (as a function object) and alpha. Note that the m(t) function must be a cdf, (must have codomain in (0,1)).


Note: Default preliminary estimates are based on standard Bass model parameters. (This starting point does not always guarantee convergence of the algorithm, and preliminary estimates may be manually selected).

}

\usage{
GG.model(sales, prelimestimates = NULL, mt = "base", alpha = 0.05, ous=100,
display=T,max.iter=100,...)
}

\arguments{
\item{sales}{Instantaneous sales}
\item{prelimestimates}{Vector containing the preliminary estimates of the parameters, default values are based on parameter estimates of a standard Bass model}
\item{mt}{Function type object, representing the variable market potential, the default for m(t) is that of GGM}
\item{alpha}{Desired significance level, the default value is 0.05}
\item{ous}{Numeric value for the out-of-sample forecasts, the default value is 100}
\item{display}{T or F, to display the plot of the model or not, the default value is T}
\item{max.iter}{Maximum number of iterations, the default value is 100}
\item{...}{Other graphic parameters}


}
\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}

\examples{
# Example 1


# dati <- c(169,397,1496,2131,2678,3431,3852,4725,5081,4592,
# 6272,6572,6479,7092,6669,7498,7380,5993,5882,9523,9885,9437
# ,10023,10103,9534,11228,10779,10687,11732,11460,12142,11465,
# 11854,11177,11112,11324,12790,12229,12116,11280,14460,13090,
# 12383,13076,13518,13781,13455,13758,14747,12405,8145,11245,
# 12211,14557,13943,14838,14275,14911,14003,14111,14241,13242,
# 15477,15219,14691,14541,12465,15909,16118,10568,11235,17345,
# 15694,15746,17129,16127,15691,16689,16552,16326,16485,15615,
# 17040,16119,13731,16102,14692,14162,17013,17058,15782,14762,
# 16813,16152,15954,16129,16356,16752)



# sp = c(1.69062e+06,2.60513e-03,3.20522e-02,1.00000e-03,1.00000e-01)
# sp1 = c(1.69062e+06,2.60513e-03,3.20522e-02)

# GG.model(sales = dati, prelimestimates = sp1 , function(x) pchisq(x,10),col=2)
# GG.model(sales = dati, mt = function(x) pchisq(x,10),col=2)
# GG.model(sales = dati, prelimestimates = sp,col=2)
# GG.model(sales = dati, col=2)
}

\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{BASS.plot}}

\code{\link{make.instantaneous}}

\code{\link{BASS.generalized}}

\code{\link{UCRCD}}


}

\keyword{GUSEO GUIDLIN}
