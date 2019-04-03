\name{LogisticLocationScaleFamily}
\alias{LogisticLocationScaleFamily}
\alias{LOGISTINT2}


\title{Generating function for Logistic location and scale families}
\description{
  Generates an object of class \code{"L2LocationScaleFamily"} which
  represents a normal location and scale family.
}
\usage{
LogisticLocationScaleFamily(location = 0, scale = 1, trafo)
LOGISTINT2
}
\arguments{
  \item{location}{ location }
  \item{scale}{ scale }
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
  \code{LOGISTINT2} is a constant used for the scale part of the Fisher information.
  More precisely \code{LOGISTINT2} equals to
  \eqn{\int_{-\infty}^{\infty} (\tanh(x/2)\,x-1)^2\,{\rm dlogis}(x)\,dx}{%
  integral_{-Inf}^{Inf} (tanh(x/2)x-1)^2 dlogis(x) dx}.
}
\value{Object of class \code{"L2LocationScaleFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Logis-class}}}
\examples{
(L1 <- LogisticLocationScaleFamily())
plot(L1)
FisherInfo(L1)
### need smaller integration range:
distrExoptions("ElowerTruncQuantile"=1e-4,"EupperTruncQuantile"=1e-4)
checkL2deriv(L1)
distrExoptions("ElowerTruncQuantile"=1e-7,"EupperTruncQuantile"=1e-7)
##
set.seed(123)
x <- rlogis(100,location=1,scale=2)
CvMMDEstimator(x, L1)
}
\concept{Logistic location and scale model}
\concept{location  and scale model}
\keyword{models}
