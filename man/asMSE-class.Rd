\name{asMSE-class}
\docType{class}
\alias{asMSE-class}

\title{Asymptotic mean square error}
\description{Class of asymptotic mean square error.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asMSE", ...)}.
  More frequently they are created via the generating function 
  \code{asMSE}.   
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}: 
      \dQuote{asymptotic mean square error}. }
    \item{\code{biastype}:}{Object of class \code{"BiasType"}: 
      symmetric, one-sided or asymmetric }
    \item{\code{normtype}:}{Object of class \code{"NormType"}: 
      norm in which a multivariate parameter is considered}
  }
}
\section{Extends}{
Class \code{"asGRisk"}, directly.\cr
Class \code{"asRiskwithBias"}, by class \code{"asGRisk"}.\cr
Class \code{"asRisk"}, by class \code{"asRiskwithBias"}.\cr
Class \code{"RiskType"}, by class \code{"asGRisk"}.
}
\section{Methods}{
No methods defined with class "asMSE" in the signature.
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asGRisk-class}}, \code{\link{asMSE}}}
\examples{
new("asMSE")
}
\concept{mean square error}
\concept{risk}
\keyword{classes}
