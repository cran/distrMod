\name{NormScaleFamily}
\alias{NormScaleFamily}

\title{Generating function for normal scale families}
\description{
  Generates an object of class \code{"L2ScaleFamily"} which
  represents a normal scale family.
}
\usage{
NormScaleFamily(sd = 1, mean = 0, trafo)
}
\arguments{
  \item{sd}{ standard deviation }
  \item{mean}{ mean }
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
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
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Norm-class}}}
\examples{
(N1 <- NormScaleFamily())
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)
}
\concept{normal scale model}
\concept{scale model}
\keyword{models}
