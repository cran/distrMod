\name{NormLocationFamily}
\alias{NormLocationFamily}
\title{Generating function for normal location families}
\description{
  Generates an object of class \code{"L2LocationFamily"} which
  represents a normal location family.
}
\usage{
NormLocationFamily(mean = 0, sd = 1, trafo)
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
\value{Object of class \code{"L2LocationFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Norm-class}}}
\examples{
(N1 <- NormLocationFamily())
plot(N1)
L2derivDistr(N1)
}
\concept{normal location model}
\concept{location model}
\keyword{models}
