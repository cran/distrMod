\name{NormType}
\alias{NormType}

\title{Generating function for NormType-class}
\description{
  Generates an object of class \code{"NormType"}.
}
\usage{NormType(name = "EuclideanNorm", fct = EuclideanNorm)}
\arguments{
  \item{name}{ slot \code{name} of the class}
  \item{fct}{  slot \code{fct} of the class }
}
%\details{}
\value{Object of class \code{"NormType"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{NormType-class}}}
\examples{
NormType()
}
\concept{asymptotic covariance}
\concept{risk}
\keyword{math}
