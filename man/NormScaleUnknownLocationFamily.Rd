\name{NormScaleUnknownLocationFamily}
\alias{NormScaleUnknownLocationFamily}

\title{Generating function for normal scale families with unknown location as nuisance}
\description{
  Generates an object of class \code{"L2LocationScaleFamily"} which
  represents a normal scale family with unknown location as nuisance.
}
\usage{
NormScaleUnknownLocationFamily(sd = 1, mean = 0, trafo)
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
(N1 <- NormScaleUnknownLocationFamily())
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)
}
\concept{normal location model}
\concept{location model}
\keyword{models}
