\name{asRiskwithBias-class}
\docType{class}
\alias{asRiskwithBias-class}
\alias{biastype,asRiskwithBias-method}
\alias{biastype<-,asRiskwithBias-method}
\alias{normtype,asRiskwithBias-method}
\alias{normtype<-,asRiskwithBias-method}
\alias{norm,asRiskwithBias-method}
\alias{normtype}
\alias{normtype<-}
\alias{norm}
\alias{biastype}
\alias{biastype<-}

\title{Aymptotic risk}
\description{Class of asymptotic risks.}
\section{Objects from the Class}{A ``virtual'' Class (although it does not contain "VIRTUAL"): No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}.}
    \item{\code{biastype}:}{Object of class \code{"BiasType"}.}
    \item{\code{normtype}:}{Object of class \code{"NormType"}.}
  }
}
\section{Extends}{
Class \code{"RiskType"}, directly.
}
\section{Methods}{
  \describe{
    \item{biastype}{\code{signature(object = "asRiskwithBias")}: 
      accessor function for slot \code{biastype}. }
    \item{biastype<-}{\code{signature(object = "asRiskwithBias", value = "BiasType")}: 
      replacement function for slot \code{biastype}. }
      accessor function for slot \code{fct} of the slot \code{normtype}. }
    \item{normtype}{\code{signature(object = "asRiskwithBias")}: 
      accessor function for slot \code{normtype}. }
    \item{normtype<-}{\code{signature(object = "asRiskwithBias", value = "NormType")}: 
      replacement function for slot \code{normtype}. }
    \item{norm}{\code{signature(object = "asRiskwithBias")}: 
        }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{asRisk-class}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}
