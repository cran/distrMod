\name{CauchyLocationFamily}
\alias{CauchyLocationFamily}

\title{Generating function for Cauchy location families}
\description{
  Generates an object of class \code{"L2LocationFamily"} which
  represents a Cauchy location family.
}
\usage{
CauchyLocationFamily(loc = 0, scale = 1, trafo)
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
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Cauchy-class}}}
\examples{
(C1 <- CauchyLocationFamily())
plot(C1)
FisherInfo(C1)
### need smaller integration range:
checkL2deriv(C1)
}
\concept{Cauchy location and scale model}
\concept{location  and scale model}
\keyword{models}
