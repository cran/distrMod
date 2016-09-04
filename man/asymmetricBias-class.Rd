\name{asymmetricBias-class}
\docType{class}
\alias{asymmetricBias-class}
\alias{nu,asymmetricBias-method}
\alias{nu<-,asymmetricBias-method}
\alias{nu}
\alias{nu<-}

\title{asymmetric Bias Type}
\description{Class of asymmetric bias types.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asymmetricBias", ...)}.
  More frequently they are created via the generating function 
  \code{asymmetricBias}.   
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}.}
    \item{\code{nu}}{Object of class \code{"numeric"}; 
                      to be in (0,1] x (0,1] with maximum 1; weights for 
                      negative and positive bias, respectively}
  }
}
\section{Methods}{
  \describe{
    \item{nu}{\code{signature(object = "asymmetricBias")}: 
      accessor function for slot \code{nu}. }
    \item{nu<-}{\code{signature(object = "asymmetricBias", value = "numeric")}: 
      replacement function for slot \code{nu}. }
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
asymmetricBias()
## The function is currently defined as
function(){ new("asymmetricBias", name = "asymmetric Bias", nu = c(1,1)) }

aB <- asymmetricBias()
nu(aB)
try(nu(aB) <- -2) ## error
nu(aB) <- c(0.3,1)
}
\concept{asymptotic risk}
\concept{bias}
\concept{risk}
\keyword{classes}
