\name{asCov}
\alias{asCov}

\title{Generating function for asCov-class}
\description{
  Generates an object of class \code{"asCov"}.
}
\usage{asCov()}
%\details{}
\value{Object of class \code{"asCov"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asCov-class}}}
\examples{
asCov()

## The function is currently defined as
function(){ new("asCov") }
}
\concept{asymptotic covariance}
\concept{risk}
\keyword{robust}
