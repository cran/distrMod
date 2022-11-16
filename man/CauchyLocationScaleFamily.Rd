\name{CauchyLocationScaleFamily}
\alias{CauchyLocationScaleFamily}
\alias{CauchyFamily}

\title{Generating function for Cauchy location and scale families}
\description{
  Generates an object of class \code{"L2LocationScaleFamily"} which
  represents a Cauchy location and scale family.
}
\usage{
CauchyLocationScaleFamily(loc = 0, scale = 1, trafo)
}
\arguments{
  \item{loc}{ location }
  \item{scale}{ scale }
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"L2LocationScaleFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Cauchy-class}}}
\examples{
(C1 <- CauchyLocationScaleFamily())
## synonymous: C1 <- CauchyFamily()
plot(C1)
FisherInfo(C1)
### need smaller integration range:
distrExoptions("ElowerTruncQuantile"=1e-4,"EupperTruncQuantile"=1e-4)
checkL2deriv(C1)
distrExoptions("ElowerTruncQuantile"=1e-7,"EupperTruncQuantile"=1e-7)
}
\concept{Cauchy location and scale model}
\concept{location  and scale model}
\keyword{models}
