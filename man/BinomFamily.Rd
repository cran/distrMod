\name{BinomFamily}
\alias{BinomFamily}

\title{Generating function for Binomial families}
\description{
  Generates an object of class \code{"L2ParamFamily"} which
  represents a Binomial family where the probability of
  success is the parameter of interest.
}
\usage{
BinomFamily(size = 1, prob = 0.5, trafo)
}
\arguments{
  \item{size}{ number of trials }
  \item{prob}{ probability of success }
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"L2ParamFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the 
  Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Binom-class}}}
\examples{
(B1 <- BinomFamily(size = 25, prob = 0.25))
plot(B1)
FisherInfo(B1)
checkL2deriv(B1)
}
\concept{binomial model}
\keyword{models}
