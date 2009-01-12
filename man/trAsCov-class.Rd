\name{trAsCov-class}
\docType{class}
\alias{trAsCov-class}
\title{Trace of asymptotic covariance}
\description{Class of trace of asymptotic covariance.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("trAsCov", ...)}.
  More frequently they are created via the generating function 
  \code{trAsCov}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{trace of asymptotic covariance}. }
  }
}
\section{Extends}{
Class \code{"asRisk"}, directly.\cr
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
\seealso{\code{\link{asRisk-class}}, \code{\link{trAsCov}}}
\examples{
new("trAsCov")
}
\concept{asymptotic covariance}
\concept{risk}
\keyword{classes}
