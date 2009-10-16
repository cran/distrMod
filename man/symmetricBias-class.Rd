\name{symmetricBias-class}
\docType{class}
\alias{symmetricBias-class}

\title{symmetric Bias Type}
\description{Class of symmetric bias types.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("symmetricBias", ...)}.
  More frequently they are created via the generating function 
  \code{symmetricBias}.   
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}.}
  }
}
\section{Methods}{
No methods defined with class "symmetricBias" in the signature.
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
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{BiasType-class}}}
\examples{
symmetricBias()
## The function is currently defined as
function(){ new("symmetricBias", name = "symmetric Bias") }
}
\concept{asymptotic risk}
\concept{bias}
\concept{risk}
\keyword{classes}
