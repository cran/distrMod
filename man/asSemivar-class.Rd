\name{asSemivar-class}
\docType{class}
\alias{asSemivar-class}
\alias{sign,asSemivar-method}
\alias{sign<-,asSemivar-method}

\title{Semivariance Risk Type}
\description{Class for semi-variance risk.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asSemivar", ...)}.
  More frequently they are created via the generating function 
  \code{asSemivar}.   
}
\section{Slots}{
  \describe{
    \item{\code{type}}{Object of class \code{"character"}: 
      \dQuote{asymptotic mean square error}. }
    \item{\code{biastype}}{Object of class \code{"BiasType"}: 
      symmetric, one-sided or asymmetric }
    \item{\code{normtype}}{Object of class \code{"NormType"}: 
      norm in which a multivariate parameter is considered}
  }
}
\section{Methods}{
  \describe{
    \item{sign}{\code{signature(object = "asSemivar")}:
      accessor function for slot \code{sign}. }
    \item{sign<-}{\code{signature(object = "asSemivar", value = "numeric")}:
      replacement function for slot \code{sign}. }
        }
}

\section{Extends}{
Class \code{"asGRisk"}, directly.\cr
Class \code{"asRiskwithBias"}, by class \code{"asGRisk"}.\cr
Class \code{"asRisk"}, by class \code{"asRiskwithBias"}.\cr
Class \code{"RiskType"}, by class \code{"asGRisk"}.
}
\references{
  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{asGRisk-class}}, \code{\link{asMSE}}}
\examples{
asSemivar()
}
\concept{asymptotic risk}
\concept{semivariance}
\concept{bias}
\concept{risk}
\keyword{classes}
