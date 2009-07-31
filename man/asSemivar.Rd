\name{asSemivar}
\alias{asSemivar}

\title{Generating function for asSemivar-class}
\description{
  Generates an object of class \code{"asSemivar"}.
}
\usage{asSemivar(sign = 1)}
%\details{}
\arguments{
  \item{sign}{ positive (=1) or negative Bias (=-1)}
}

\value{Object of class \code{"asSemivar"}}
\references{
  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{onesidedBias-class}}}
\examples{
asSemivar()

}
\concept{bias}
\concept{semivariance}
\keyword{robust}
