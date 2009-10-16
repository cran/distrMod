\name{BiasType-class}
\docType{class}
\alias{BiasType-class}
\alias{name,BiasType-method}
\alias{name<-,BiasType-method}

\title{Bias Type}
\description{Class of bias types.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}.}
  }
}
\section{Methods}{
  \describe{
    \item{name}{\code{signature(object = "BiasType")}: 
      accessor function for slot \code{name}. }
    \item{name<-}{\code{signature(object = "BiasType", value = "character")}: 
      replacement function for slot \code{name}. }
        }
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
\seealso{\code{\link{RiskType-class}}}
\examples{
aB <- positiveBias()
name(aB)
}
\concept{asymptotic risk}
\concept{bias}
\concept{risk}
\keyword{classes}
