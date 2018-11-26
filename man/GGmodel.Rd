\name{GG.model}
\alias{GG.model}
\title{Function that returns the results of Guseo-Guidolin model}
\description{Function returns the results of dynamic market potential model, with m(t) defined as the classic Guseo-Guidolin model or defining yout own m(t) function (which must depend only on t, and must be similiar to a repartition function).
To use this function, you have two choices:

1) Use the classic GG model, and you have to pass to the function: sales, preliminar estimates and alpha

2) Use the model definig your own m(t) function, in this case you have to pass to the function: sales,preliminar estimates, mt (as a function object) and alpha.
Remember that the m(t) function must be positive and similar to a repartition function, in fact, the function must have codomain in (0,1).

N.B.

You don't have to pass the preliminar estimates because the default values are taken by running the BASS.standard function.
This doesn't guarantee that with this starting point the algorith converge!
}

\usage{
GG.model(sales,prelimestimates,mt,alpha=0.05,...)
}

\arguments{
\item{sales}{The istantanuous sales}
\item{prelimestimates}{A vector containing the preliminar estimates of the parameters, the default value is taken by BASS standard estimate}
\item{mt}{A function type object, which represent the variable market, if you don't specify a function the default is the classic GG-model}
\item{alpha}{The desired significativity for the interval of confidence, the defaul value is 0.05}
\item{...}{Other grafics parameters}

}
\author{
Zanghi Federico \email{federico.zanghi@studenti.unipd.it}
Department of Statistics, University of Padova, Italy \url{http://www.unipd.it}
}

\examples{
# Example 1


dati <- c(169,397,1496,2131,2678,3431,3852,4725,5081,4592,
6272,6572,6479,7092,6669,7498,7380,5993,5882,9523,9885,9437
,10023,10103,9534,11228,10779,10687,11732,11460,12142,11465,
11854,11177,11112,11324,12790,12229,12116,11280,14460,13090,
12383,13076,13518,13781,13455,13758,14747,12405,8145,11245,
12211,14557,13943,14838,14275,14911,14003,14111,14241,13242,
15477,15219,14691,14541,12465,15909,16118,10568,11235,17345,
15694,15746,17129,16127,15691,16689,16552,16326,16485,15615,
17040,16119,13731,16102,14692,14162,17013,17058,15782,14762,
16813,16152,15954,16129,16356,16752)



sp = c(1.69062e+06,2.60513e-03,3.20522e-02,1.00000e-03,1.00000e-01)
sp1 = c(1.69062e+06,2.60513e-03,3.20522e-02)

GG.model(sales = dati, prelimestimates = sp1 , function(x) pchisq(x,10),col=2)
GG.model(sales = dati, mt = function(x) pchisq(x,10),col=2)
GG.model(sales = dati, prelimestimates = sp,col=2)
GG.model(sales = dati, col=2)
}

\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{BASS.plot}}

\code{\link{make.instantaneous}}

\code{\link{BASS.generalized}}


}

\keyword{GUSEO GUIDLIN}
