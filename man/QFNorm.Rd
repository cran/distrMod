\name{QFNorm}
\alias{QFNorm}

\title{Generating function for QFNorm-class}
\description{
  Generates an object of class \code{"QFNorm"}.
}
\usage{QFNorm(name = "norm based on quadratic form", 
       QuadForm = PosSemDefSymmMatrix(matrix(1)))}

\arguments{
  \item{name}{ slot \code{name} of the class}
  \item{QuadForm}{  slot \code{QuadForm} of the class }
}
\value{Object of class \code{"QFNorm"}}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{QFNorm-class}}}
\examples{
## IGNORE_RDIFF_BEGIN
QFNorm()

## The function is currently defined as
function(){ new("QFNorm") }
## IGNORE_RDIFF_END
}
\concept{asymptotic covariance}
\concept{risk}
\keyword{math}
