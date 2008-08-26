\name{ExpScaleFamily}
\alias{ExpScaleFamily}

\title{Generating function for exponential scale families}
\description{
  Generates an object of class \code{"L2ScaleFamily"} which
  represents an exponential scale family.
}
\usage{
ExpScaleFamily(scale = 1, trafo)
}
\arguments{
  \item{scale}{ scale (= 1/rate) }
  \item{trafo}{ function in \code{param} or matrix: optional transformation of the parameter }
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled. The scale parameter corresponds
  to \eqn{1/\code{rate}}.
}
\value{Object of class \code{"L2ScaleFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Exp-class}}}
\examples{
(E1 <- ExpScaleFamily())
plot(E1)
Map(L2deriv(E1)[[1]])
checkL2deriv(E1)
}
\concept{exponential scale model}
\concept{scale model}
\keyword{models}
