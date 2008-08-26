\name{negativeBias}
\alias{negativeBias}

\title{Generating function for onesidedBias-class}
\description{
  Generates an object of class \code{"onesidedBias"}.
}
\usage{negativeBias(name = "negative Bias")}
%\details{}
\arguments{
  \item{name}{ name of the bias type}
}

\value{Object of class \code{"onesidedBias"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{onesidedBias-class}}}
\examples{
negativeBias()

## The function is currently defined as
function(){ new("onesidedBias", name = "negative Bias", sign = -1) }
}
\concept{bias}
\concept{risk}
\keyword{robust}
