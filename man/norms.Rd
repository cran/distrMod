\name{norm}
\alias{EuclideanNorm}
\alias{QuadFormNorm}

\title{Norm functions}
\description{
  Functions to determine certain norms.
}
\usage{EuclideanNorm(x)
QuadFormNorm(x,A)}
%\details{}
\arguments{
  \item{x}{vector or matrix; norm is determined columnwise}
  \item{A}{pos. semidefinite Matrix}
}

\value{the columnwise evaluated norms}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{onesidedBias-class}}}
\examples{
mm <- matrix(rnorm(20),2,10)
EuclideanNorm(mm)
QuadFormNorm(mm, A = PosSemDefSymmMatrix(matrix(c(3,1,1,1),2,2)))
}
\concept{bias}
\concept{risk}
\keyword{robust}
