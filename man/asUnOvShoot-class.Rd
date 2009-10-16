\name{asUnOvShoot-class}
\docType{class}
\alias{asUnOvShoot-class}
\alias{width}
\alias{width,asUnOvShoot-method}
\alias{show,asUnOvShoot-method}

\title{Asymptotic under-/overshoot probability}
\description{Class of asymptotic under-/overshoot probability.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asUnOvShoot", ...)}.
  More frequently they are created via the generating function 
  \code{asUnOvShoot}.   
}
\section{Slots}{
  \describe{
    \item{\code{type}}{Object of class \code{"character"}:
      \dQuote{asymptotic under-/overshoot probability}. }
    \item{\code{width}}{Object of class \code{"numeric"}: 
      half the width of given confidence interval. }
    \item{\code{biastype}}{Object of class \code{"BiasType"}: 
      symmetric, one-sided or asymmetric }
  }
}
\section{Extends}{
Class \code{"asGRisk"}, directly.\cr
Class \code{"asRiskwithBias"}, by class \code{"asGRisk"}.\cr
Class \code{"asRisk"}, by class \code{"asRiskwithBias"}.\cr
Class \code{"RiskType"}, by class \code{"asGRisk"}.
}
\section{Methods}{
  \describe{
    \item{width}{\code{signature(object = "asUnOvShoot")}: 
      accessor function for slot \code{width}. }
    \item{show}{\code{signature(object = "asUnOvShoot")}}
  }
}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asGRisk-class}}}
\examples{
new("asUnOvShoot")
}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}
