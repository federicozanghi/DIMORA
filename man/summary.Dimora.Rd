\name{summary.Dimora}
\alias{summary.Dimora}
\title{summary method for Dimora objects

}
\description{
  Summary method for the \code{Dimora} objects.
}
\usage{
\method{summary}{Dimora}(object,...)
}
\arguments{
  \item{object}{An object of class \code{Dimora} i.e.  \code{BASS.generalized}, \code{BASS.standard}, \code{GG.model} and \code{UCRCD}.}
  \item{...}{Additional arguments related to the summary produced.}
}

\value{
  \code{summary} is a generic function used to produce summary of results for the models fitted in Dimora. The function invokes particular methods which depend on the class of the first argument, in this case the class is \code{Dimora}.}


\author{
Zanghi Federico \email{federico.zanghi.11@gmail.com}}
