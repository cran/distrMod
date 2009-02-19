\name{fiCov}
\alias{fiCov}

\title{Generating function for fiCov-class}
\description{
  Generates an object of class \code{"fiCov"}.
}
\usage{asCov()}
%\details{}
\value{Object of class \code{"fiCov"}}
\references{
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiCov-class}}}
\examples{
fiCov()

## The function is currently defined as
function(){ new("fiCov") }
}
\concept{finite-sample covariance}
\concept{risk}
\keyword{robust}
