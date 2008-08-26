\name{trFiCov-class}
\docType{class}
\alias{trFiCov-class}

\title{Trace of finite-sample covariance}
\description{Class of trace of finite-sample covariance.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("trFiCov", ...)}.
  More frequently they are created via the generating function 
  \code{trFiCov}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{trace of finite-sample covariance}. }
  }
}
\section{Extends}{
Class \code{"fiRisk"}, directly.\cr
Class \code{"RiskType"}, by class \code{"fiRisk"}.
}
%\section{Methods}{}
\references{
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiRisk-class}}, \code{\link{trFiCov}}}
\examples{
new("trFiCov")
}
\concept{finite-sample covariance}
\concept{risk}
\keyword{classes}
