\name{NormLocationUnknownScaleFamily}
\alias{NormLocationUnknownScaleFamily}

\title{Generating function for normal location families with unknown scale as nuisance}
\description{
  Generates an object of class \code{"L2LocationScaleFamily"} which
  represents a normal location family with unknown scale as nuisance.
}
\usage{
NormLocationUnknownScaleFamily(mean = 0, sd = 1, trafo)
}
\arguments{
  \item{mean}{ mean }
  \item{sd}{ standard deviation }
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
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Norm-class}}}
\examples{
(N1 <- NormLocationUnknownScaleFamily())
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)
}
\concept{normal location model}
\concept{location model}
\keyword{models}
