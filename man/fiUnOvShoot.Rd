\name{fiUnOvShoot}
\alias{fiUnOvShoot}

\title{Generating function for fiUnOvShoot-class}
\description{
  Generates an object of class \code{"fiUnOvShoot"}.
}
\usage{fiUnOvShoot(width = 1.960)}
\arguments{
  \item{width}{ positive real: half the width of given confidence interval. }
}
%\details{}
\value{Object of class \code{"fiUnOvShoot"}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.
  
  Rieder, H. (1989) A finite-sample minimax regression estimator. 
  Statistics \bold{20}(2): 211--221.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Ruckdeschel, P. and Kohl, M. (2005) How to approximate 
  the finite sample risk of M-estimators.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiUnOvShoot-class}}}
\examples{
fiUnOvShoot()

## The function is currently defined as
function(width = 1.960){ new("fiUnOvShoot", width = width) }
}
\concept{finite-sample risk}
\concept{risk}
\keyword{robust}
