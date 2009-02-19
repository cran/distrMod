\name{fiUnOvShoot-class}
\docType{class}
\alias{fiUnOvShoot-class}
\alias{width,fiUnOvShoot-method}
\alias{show,fiUnOvShoot-method}

\title{Finite-sample under-/overshoot probability}
\description{Class of finite-sample under-/overshoot probability.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("fiUnOvShoot", ...)}.
  More frequently they are created via the generating function 
  \code{fiUnOvShoot}.   
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      dQuote{finite-sample under-/overshoot probability}. }
    \item{\code{width}:}{Object of class \code{"numeric"}: 
      half the width of given confidence interval. }
  }
}
\section{Extends}{
Class \code{"fiRisk"}, directly.\cr
Class \code{"RiskType"}, by class \code{"fiRisk"}.
}
\section{Methods}{
  \describe{
    \item{width}{\code{signature(object = "fiUnOvShoot")}: 
      accessor function for slot \code{width}. }
    \item{show}{\code{signature(object = "fiUnOvShoot")}}
  }
}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.
  
  Rieder, H. (1989) A finite-sample minimax regression estimator. 
  Statistics \bold{20}(2): 211--221.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Ruckdeschel, P. and Kohl, M. (2005) Computation of the Finite Sample Risk 
  of M-estimators on Neighborhoods.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiRisk-class}}}
\examples{
new("fiUnOvShoot")
}
\concept{finite-sample risk}
\concept{risk}
\keyword{classes}
