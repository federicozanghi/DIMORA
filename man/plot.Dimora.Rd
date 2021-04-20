\name{plot.Dimora}
\alias{plot.Dimora}
\title{plot method for Dimora}
\description{
 Three plots (to be selected with \code{typeog}) are currently available: a plot of fitted values over observed data (typeog=1), a plot of residuals against fitted values (typeog=2), a plot of residuals against time (typeog=3).

}
\usage{
\method{plot}{Dimora}(x,..., product=1, mode=c("i","c"),typeog=1,oos=c("y","n"),
                    xlim=NULL,ltycurve=1,lwdcurve=1.5,typecurve="l",colcurve=2)
}
\arguments{
  \item{x}{An object of class \code{Dimora}, i.e.  \code{BASS.generalized}, \code{BASS.standard}, \code{GG.model}.}
  \item{product}{1 or 2, which indicates whether there are one or two products}
  \item{mode}{Mode "i" for instantaneous values and "c" for cumulative values}
  \item{typeog}{Fitted values on observed data \code{typeog=1}, residuals against fitted values \code{typeog=2} and residuals against time \code{typeog=3}, and all the plots with \code{typeog=0}}
  \item{oos}{"y" or "n" if the out-of-sample plot is desired}
  \item{xlim}{The limit for the graph, default values is null, if the out-of-sample plot is desired the portion of out-of-sample values may be defined}
  \item{ltycurve}{"lty" parameter for the fitted values plot (\code{typeog=1}).}
  \item{lwdcurve}{"lwd" parameter for the fitted values plot (\code{typeog=1}).}
  \item{typecurve}{"type" parameter for the fitted values plot (\code{typeog=1}).}
  \item{colcurve}{"col" parameter for the fitted values plot (\code{typeog=1}).}
  \item{...}{Other graphic parameters.}
}

\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}
