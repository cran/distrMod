\name{symmetricBias}
\alias{symmetricBias}

\title{Generating function for symmetricBias-class}
\description{
  Generates an object of class \code{"symmetricBias"}.
}
\usage{symmetricBias(name = "symmetric Bias")}
%\details{}
\arguments{
  \item{name}{ name of the bias type}
}

\value{Object of class \code{"symmetricBias"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{symmetricBias-class}}}
\examples{
symmetricBias()

## The function is currently defined as
function(){ new("symmetricBias", name = "symmetric Bias") }
}
\concept{bias}
\concept{risk}
\keyword{robust}
