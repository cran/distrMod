\name{asHampel}
\alias{asHampel}

\title{Generating function for asHampel-class}
\description{
  Generates an object of class \code{"asHampel"}.
}
\usage{asHampel(bound = Inf, biastype = symmetricBias(), normtype = NormType())}
\arguments{
  \item{bound}{ positive real: bias bound }
  \item{biastype}{ a bias type of class \code{BiasType}}
  \item{normtype}{ a norm type of class \code{NormType}}
}
%\details{}
\value{Object of class \code{asHampel}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asHampel-class}}}
\examples{
asHampel()

## The function is currently defined as
function(bound = Inf, biastype = symmetricBias(), normtype = NormType()){ 
    new("asHampel", bound = bound, biastype = biastype, normtype = normtype) }
}
\concept{Hampel risk}
\concept{risk}
\keyword{robust}
