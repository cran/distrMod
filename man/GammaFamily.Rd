\name{GammaFamily}
\alias{GammaFamily}

\title{Generating function for Gamma families}
\description{
  Generates an object of class \code{"L2ParamFamily"} which
  represents a Gamma family.
}
\usage{
GammaFamily(scale = 1, shape = 1, trafo, withL2derivDistr = TRUE)
}
\arguments{
  \item{scale}{ positive real: scale parameter }
  \item{shape}{ positive real: shape parameter }
  \item{trafo}{ matrix: transformation of the parameter }
  \item{withL2derivDistr}{logical: shall the distribution of the L2 derivative
      be computed? Defaults to \code{TRUE}; setting it to \code{FALSE}
      speeds up computations.}
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"L2ParamFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Gammad-class}}}
\examples{
(G1 <- GammaFamily())
FisherInfo(G1)
checkL2deriv(G1)
}
\concept{Gamma model}
\keyword{models}
