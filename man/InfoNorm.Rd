\name{InfoNorm}
\alias{InfoNorm}

\title{Generating function for InfoNorm-class}
\description{
  Generates an object of class \code{"InfoNorm"} ---
  used for information-standardized influence curves.
}
\usage{InfoNorm()}
%\details{}
\value{Object of class \code{"InfoNorm"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfoNorm-class}}}
\examples{
## IGNORE_RDIFF_BEGIN
InfoNorm()

## The function is currently defined as
function(){ new("InfoNorm") }
## IGNORE_RDIFF_END
}
\concept{asymptotic covariance}
\concept{risk}
\keyword{robust}
