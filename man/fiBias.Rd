\name{fiBias}
\alias{fiBias}

\title{Generating function for fiBias-class}
\description{
  Generates an object of class \code{"fiBias"}.
}
\usage{fiBias()}
%\details{}
\value{Object of class \code{"fiBias"}}
\references{
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiBias-class}}}
\examples{
fiBias()

## The function is currently defined as
function(){ new("fiBias") }
}
\concept{finite-sample bias}
\concept{risk}
\keyword{robust}
