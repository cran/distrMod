\name{PoisFamily}
\alias{PoisFamily}

\title{Generating function for Poisson families}
\description{
  Generates an object of class \code{"L2ParamFamily"} which
  represents a Poisson family.
}
\usage{
PoisFamily(lambda = 1, trafo)
}
\arguments{
  \item{lambda}{ positive mean }
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"L2ParamFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Pois-class}}}
\examples{
(P1 <- PoisFamily(lambda = 4.5))
plot(P1)
FisherInfo(P1)
checkL2deriv(P1)
}
\concept{Poisson model}
\keyword{models}
