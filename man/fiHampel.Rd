\name{fiHampel}
\alias{fiHampel}

\title{Generating function for fiHampel-class}
\description{
  Generates an object of class \code{"fiHampel"}.
}
\usage{fiHampel(bound = Inf)}
\arguments{
  \item{bound}{ positive real: bias bound }
}
%\details{}
\value{Object of class \code{fiHampel}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators.   
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiHampel-class}}}
\examples{
fiHampel()

## The function is currently defined as
function(bound = Inf){ new("fiHampel", bound = bound) }
}
\concept{Hampel risk}
\concept{risk}
\keyword{robust}
