\name{UCRCD}
\alias{UCRCD}
\title{Function that estimates a UCRCD model}
\description{This function allows to estimate diachronic competition between two products, according to the Unbalanced Competition Regime Change Diachronic model, UCRCD, (Guseo and Mortarino, 2014).}

\usage{
UCRCD(sales1, sales2, c2, display=T, alpha=0.05,
      delta=0.01, gamma=0.01, par="double",
      m1 =BASS.standard(sales1,display = F)$Estimate[1,1],
      m2 =BASS.standard(sales2,display = F)$Estimate[1,1],
      p1c=BASS.standard(sales1,display = F)$Estimate[2,1],
      q1c=BASS.standard(sales1,display = F)$Estimate[3,1],
      p2c=BASS.standard(sales2,display = F)$Estimate[2,1],
      q2c=BASS.standard(sales2,display = F)$Estimate[3,1])
}

\arguments{
\item{sales1}{Instantaneous sales of the first product}
\item{sales2}{Instantaneous sales of the second product}
\item{c2}{Time at which the second product enters the market}
\item{display}{T or F, to display the plot of the model or not, the default value is T}
\item{alpha}{Desired significance level, the default value is 0.05}
\item{delta}{Preliminary estimates for parameter delta}
\item{gamma}{Preliminary estimates for parameter gamma}
\item{par}{Default value is par = ”double”, where delta is different from gamma. If par = "unique", delta is equal to gamma}
\item{m1}{Preliminary estimates of the first product's market potential under competition}
\item{m2}{Preliminary estimates of the second product's market potential}
\item{p1c}{Preliminary estimates of the first product's innovation coefficient under competition}
\item{p2c}{Preliminary estimates of the second product's innovation coefficient under competition}
\item{q1c}{Preliminary estimates of the first product's imitation coefficient under competition}
\item{q2c}{Preliminary estimates of the second product's imitation coefficient under competition}
}

\author{
Savio Andrea \email{svandr97@gmail.com}
}


\seealso{

\code{\link{BASS.standard}}

\code{\link{BASS.standard.generator}}

\code{\link{BASS.generalized}}

\code{\link{BASS.plot}}

\code{\link{make.instantaneous}}

\code{\link{GG.model}}

}

\keyword{UCRCD}
