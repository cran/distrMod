\name{fiHampel-class}
\docType{class}
\alias{fiHampel-class}
\alias{bound,fiHampel-method}
\alias{show,fiHampel-method}

\title{Finite-sample Hampel risk}
\description{Class of finite-sample Hampel risk which is
  the trace of the finite-sample covariance subject to
  a given bias bound (bound on gross error sensitivity). }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("fiHampel", ...)}.
  More frequently they are created via the generating function 
  \code{fiHampel}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}: 
      \dQuote{trace of finite-sample covariance for given bias bound}. }
    \item{\code{bound}:}{Object of class \code{"numeric"}: 
      given positive bias bound. }
  }
}
\section{Extends}{
Class \code{"fiRisk"}, directly.\cr
Class \code{"RiskType"}, by class \code{"fiRisk"}.
}
\section{Methods}{
  \describe{
    \item{bound}{\code{signature(object = "fiHampel")}: 
      accessor function for slot \code{bound}. }
    \item{show}{\code{signature(object = "fiHampel")}}
  }
}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators.   
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiRisk-class}}, \code{\link{fiHampel}}}
\examples{
new("fiHampel")
}
\concept{Hampel risk}
\concept{risk}
\keyword{classes}
