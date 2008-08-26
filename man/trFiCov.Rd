\name{trFiCov}
\alias{trFiCov}
\title{Generating function for trFiCov-class}
\description{
  Generates an object of class \code{"trFiCov"}.
}
\usage{trFiCov()}
%\details{}
\value{Object of class \code{"trFiCov"}}
\references{
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{trFiCov-class}}}
\examples{
trFiCov()

## The function is currently defined as
function(){ new("trFiCov") }
}
\concept{finite-sample covariance}
\concept{risk}
\keyword{robust}
