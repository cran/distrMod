\name{QFNorm-class}
\docType{class}
\alias{QFNorm-class}
\alias{InfoNorm-class}
\alias{SelfNorm-class}
\alias{QuadForm,QFNorm-method}
\alias{QuadForm<-,QFNorm-method}
\alias{QuadForm}
\alias{QuadForm<-}

\title{Norm classes for norms based on quadratic forms}
\description{Classes for norms based on quadratic forms}
\section{Objects from the Class}{could be created by a call to \code{new}, but 
normally one would use the generating functions \code{QFNorm}, 
\code{InfoNorm}, and \code{SelfNorm} }

\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}.}
    \item{\code{fct}}{Object of class \code{"function"}.}
    \item{\code{QuadForm}}{Object of class \code{"PosSemDefSymmMatrix"}.}
  }
}
\section{Extends}{
\code{"QFNorm"} extends class \code{"NormType"}, directly,
and \code{"InfoNorm"} and \code{"SelfNorm"} each extend 
class \code{"QFNorm"}, directly (and do not have extra slots).
}
\section{Methods}{
  \describe{
    \item{QuadForm}{\code{signature(object = "QFNorm")}: 
      accessor function for slot \code{QuadForm}. }
    \item{QuadForm<-}{\code{signature(object = "QFNorm")}: 
      replacement function for slot \code{QuadForm}. }
        }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{NormType-class}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}
