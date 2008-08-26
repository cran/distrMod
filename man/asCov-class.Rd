\name{asCov-class}
\docType{class}
\alias{asCov-class}

\title{Asymptotic covariance}
\description{Class of asymptotic covariance.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asCov", ...)}.
  More frequently they are created via the generating function 
  \code{asCov}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{asymptotic covariance}. }
  }
}
\section{Extends}{
Class \code{"asRisk"}, directly.\cr
Class \code{"RiskType"}, by class \code{"asRisk"}.
}
\section{Methods}{
No methods defined with class "asCov" in the signature.
}
\references{ 
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asRisk-class}}, \code{\link{asCov}}}
\examples{
new("asCov")
}
\concept{asymptotic covariance}
\concept{risk}
\keyword{classes}
