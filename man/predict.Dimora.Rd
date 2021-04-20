\name{predict.Dimora}
\alias{predict.Dimora}
\title{predict method for Dimora objects

}
\description{
  Prediction of test data using any model of the package.
}
\usage{
\method{predict}{Dimora}(object,..., newx)
}
\arguments{
  \item{object}{An object of class \code{Dimora}, i.e.  \code{BASS.generalized}, \code{BASS.standard}, \code{GG.model}.}
  \item{newx}{A number or a vector of numeric values, representing time t.}
  \item{...}{Additional arguments affecting the predictions produced.}
}

\value{
  The method returns a vector of predicted values.}

\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}




