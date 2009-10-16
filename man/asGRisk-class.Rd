\name{asGRisk-class}
\docType{class}
\alias{asGRisk-class}

\title{Convex asymptotic risk}
\description{Class of special convex asymptotic risks.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}}{Object of class \code{"character"}. }
    \item{\code{biastype}}{Object of class \code{"BiasType"}: 
      symmetric, one-sided or asymmetric }
    \item{\code{normtype}}{Object of class \code{"NormType"}: 
      norm in which a multivariate parameter is considered}
  }
}
\section{Extends}{
Class \code{"asRisk"}, directly.\cr
Class \code{"RiskType"}, by class \code{"asRisk"}.
}
\section{Methods}{
No methods defined with class "asGRisk" in the signature.
}
\references{
  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asRisk-class}}}
%\examples{}
\concept{convex risk}
\concept{risk}
\keyword{classes}
