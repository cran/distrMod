\name{fiMSE}
\alias{fiMSE}

\title{Generating function for fiMSE-class}
\description{
  Generates an object of class \code{"fiMSE"}.
}
\usage{fiMSE()}
%\details{}
\value{Object of class \code{"fiMSE"}}
\references{ 
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiMSE-class}}}
\examples{
fiMSE()

## The function is currently defined as
function(){ new("fiMSE") }
}
\concept{mean square error}
\concept{risk}
\keyword{robust}
