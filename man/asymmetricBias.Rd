\name{asymmetricBias}
\alias{asymmetricBias}

\title{Generating function for asymmetricBias-class}
\description{
  Generates an object of class \code{"asymmetricBias"}.
}
\usage{asymmetricBias(name = "asymmetric Bias", nu = c(1,1) )}
%\details{}
\arguments{
  \item{name}{ name of the bias type}
  \item{nu}{weights for negative and positive bias, respectively}
}

\value{Object of class \code{"asymmetricBias"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{asymmetricBias-class}}}
\examples{
asymmetricBias()

## The function is currently defined as
function(){ new("asymmetricBias", name = "asymmetric Bias", nu = c(1,1)) }
}
\concept{bias}
\concept{risk}
\keyword{robust}
