\name{LnormScaleFamily}
\alias{LnormScaleFamily}

\title{Generating function for lognormal scale families}
\description{
  Generates an object of class \code{"L2ScaleFamily"} which
  represents a lognormal scale family.
}
\usage{
LnormScaleFamily(meanlog = 0, sdlog = 1, trafo)
}
\arguments{
  \item{meanlog}{ mean of the distribution on the log scale }
  \item{sdlog}{ standard deviation of the distribution on the log scale }
  \item{trafo}{ matrix: transformation of the parameter }
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"L2ScaleFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Lnorm-class}}}
\examples{
(L1 <- LnormScaleFamily())
plot(L1)
Map(L2deriv(L1)[[1]])
checkL2deriv(L1)
}
\concept{lognormal scale model}
\concept{scale model}
\keyword{models}
