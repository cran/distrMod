\name{asHampel-class}
\docType{class}
\alias{asHampel-class}
\alias{bound}
\alias{bound,asHampel-method}
\alias{show,asHampel-method}

\title{Asymptotic Hampel risk}
\description{Class of asymptotic Hampel risk which is
  the trace of the asymptotic covariance subject to
  a given bias bound (bound on gross error sensitivity). }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asHampel", ...)}.
  More frequently they are created via the generating function 
  \code{asHampel}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}: 
      \dQuote{trace of asymptotic covariance for given bias bound}. }
    \item{\code{bound}:}{Object of class \code{"numeric"}: 
      given positive bias bound. }
    \item{\code{biastype}:}{Object of class \code{"BiasType"}: 
      symmetric, one-sided or asymmetric }
  }
}
\section{Extends}{
Class \code{"asRiskwithBias"}, directly.\cr
Class \code{"asRisk"}, by class \code{"asRiskwithBias"}.
Class \code{"RiskType"}, by class \code{"asRisk"}.
}
\section{Methods}{
  \describe{
    \item{bound}{\code{signature(object = "asHampel")}: 
      accessor function for slot \code{bound}. }
    \item{show}{\code{signature(object = "asHampel")}}
  }
}
\references{Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asRisk-class}}, \code{\link{asHampel}}}
\examples{
new("asHampel")
}
\concept{Hampel risk}
\concept{risk}
\keyword{classes}
