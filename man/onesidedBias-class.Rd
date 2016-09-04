\name{onesidedBias-class}
\docType{class}
\alias{onesidedBias-class}
\alias{sign}
\alias{sign<-}
\alias{sign,onesidedBias-method}
\alias{sign<-,onesidedBias-method}

\title{onesided Bias Type}
\description{Class of onesided bias types.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("onesidedBias", ...)}.
  More frequently they are created via the generating function 
  \code{positiveBias} or   \code{negativeBias}.   
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}.}
    \item{\code{sign}}{Object of class \code{"numeric"}; 
                      to be in \{-1,1\} --- whether bias is to be positive or negative}
  }
}
\section{Methods}{
  \describe{
    \item{sign}{\code{signature(object = "onesidedBias")}: 
      accessor function for slot \code{sign}. }
    \item{sign<-}{\code{signature(object = "onesidedBias", value = "numeric")}: 
      replacement function for slot \code{sign}. }
        }
}
\section{Extends}{
Class \code{"BiasType"}, directly.\cr
}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{BiasType-class}}}
\examples{
positiveBias()
## The function is currently defined as
function(){ new("onesidedBias", name = "positive Bias", sign = 1) }

negativeBias()
## The function is currently defined as
function(){ new("onesidedBias", name = "negative Bias", sign = -1) }

pB <- positiveBias()
sign(pB)
try(sign(pB) <- -2) ## error
sign(pB) <- -1
}
\concept{asymptotic risk}
\concept{bias}
\concept{risk}
\keyword{classes}
