\name{fiBias-class}
\docType{class}
\alias{fiBias-class}

\title{Finite-sample Bias}
\description{Class of finite-sample bias.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("fiBias", ...)}.
  More frequently they are created via the generating function 
  \code{fiBias}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{finite-sample bias}. }
  }
}
\section{Extends}{
Class \code{"fiRisk"}, directly.\cr
Class \code{"RiskType"}, by class \code{"fiRisk"}.
}
\section{Methods}{
No methods defined with class "fiBias" in the signature.
}
\references{
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiRisk-class}}, \code{\link{fiBias}}}
\examples{
new("fiBias")
}
\concept{finite-sample bias}
\concept{risk}
\keyword{classes}
