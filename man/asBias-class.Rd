\name{asBias-class}
\docType{class}
\alias{asBias-class}

\title{Standardized Asymptotic Bias}
\description{Class of standardized asymptotic bias; i.e., 
  the neighborhood radius is omitted respectively, set to \eqn{1}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asBias", ...)}.
  More frequently they are created via the generating function 
  \code{asBias}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}}{Object of class \code{"character"}:
      \dQuote{asymptotic bias}. }
    \item{\code{biastype}}{Object of class \code{"BiasType"}: 
      symmetric, one-sided or asymmetric }
    \item{\code{normtype}}{Object of class \code{"NormType"}: 
      norm in which a multivariate parameter is considered}
  }
}
\section{Extends}{
Class \code{"asRiskwithBias"}, directly.\cr
Class \code{"asRisk"}, by class \code{"asRiskwithBias"}\cr
Class \code{"RiskType"}, by class \code{"asRisk"}.
}
%\section{Methods}{}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asRisk-class}}, \code{\link{asBias}}}
\examples{
new("asBias")
}
\concept{asymptotic bias}
\concept{risk}
\keyword{classes}
