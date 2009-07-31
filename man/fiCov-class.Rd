\name{fiCov-class}
\docType{class}
\alias{fiCov-class}

\title{Finite-sample covariance}
\description{Class of finite-sample covariance.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("fiCov", ...)}.
  More frequently they are created via the generating function 
  \code{fiCov}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{finite-sample covariance}. }
  }
}
\section{Extends}{
Class \code{"fiRisk"}, directly.\cr
Class \code{"RiskType"}, by class \code{"fiRisk"}.
}
\section{Methods}{
No methods defined with class "fiCov" in the signature.
}
\references{
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiRisk-class}}, \code{\link{fiCov}}}
\examples{
new("fiCov")
}
\concept{finite-sample covariance}
\concept{risk}
\keyword{classes}
