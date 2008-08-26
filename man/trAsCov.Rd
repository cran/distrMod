\name{trAsCov}
\alias{trAsCov}
\title{Generating function for trAsCov-class}
\description{
  Generates an object of class \code{"trAsCov"}.
}
\usage{trAsCov()}
%\details{}
\value{Object of class \code{"trAsCov"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{trAsCov-class}}}
\examples{
trAsCov()

## The function is currently defined as
function(){ new("trAsCov") }
}
\concept{asymptotic covariance}
\concept{risk}
\keyword{robust}
