\name{fiMSE-class}
\docType{class}
\alias{fiMSE-class}

\title{Finite-sample mean square error}
\description{Class of asymptotic mean square error.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("fiMSE", ...)}.
  More frequently they are created via the generating function 
  \code{fiMSE}.   
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}: 
      \dQuote{finite-sample mean square error}. }
  }
}
\section{Extends}{
Class \code{"fiRisk"}, directly.\cr
Class \code{"RiskType"}, by class \code{"fiRisk"}.
}
\section{Methods}{
No methods defined with class "fiMSE" in the signature.
}
\references{
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiRisk-class}}, \code{\link{fiMSE}}}
\examples{
new("fiMSE")
}
\concept{mean square error}
\concept{risk}
\keyword{classes}
