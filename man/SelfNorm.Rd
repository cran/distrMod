\name{SelfNorm}
\alias{SelfNorm}

\title{Generating function for SelfNorm-class}
\description{
  Generates an object of class \code{"SelfNorm"} ---
  used for self-standardized influence curves.
}
\usage{SelfNorm()}
%\details{}
\value{Object of class \code{"SelfNorm"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{SelfNorm-class}}}
\examples{
## IGNORE_RDIFF_BEGIN
SelfNorm()

## The function is currently defined as
function(){ new("SelfNorm") }
## IGNORE_RDIFF_END
}
\concept{asymptotic covariance}
\concept{risk}
\keyword{robust}
