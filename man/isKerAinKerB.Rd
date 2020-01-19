\name{isKerAinKerB}
\alias{isKerAinKerB}

\title{isKerAinKerB}
\description{
  For two matrices \code{A} and \code{B} checks whether the null space of \code{A} is a subspace
  of the null space of \code{B}, in other words, if \eqn{Ax=0} entails \code{Bx=0}.
}
\usage{
isKerAinKerB(A, B, tol = .Machine$double.eps)
}
\arguments{
  \item{A}{ a matrix; if  \code{A} is a vector, \code{A} is coerced to a matrix by \code{as.matrix}. }
  \item{B}{ a matrix; if  \code{B} is a vector, \code{B} is coerced to a matrix by \code{as.matrix}. }
  \item{tol}{ the tolerance for detecting linear dependencies in the
          columns of \code{a} and up to which the two projectors are seen as equal (see below). }
}
\details{
  via calls to \code{\link[base]{svd}}, the projectors \eqn{\pi_A}{pi_A} and \eqn{\pi_B}{pi_B}
  onto the respective orthogonal complements of \eqn{{\rm ker}(A)}{ker(A)} and
  \eqn{{\rm ker}(B)}{ker(B)} are calculated and then is checked whether 
  \eqn{\pi_B\pi_A=\pi_B}{pi_B pi_A = pi_B}.
}
\value{logical}
%\references{}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
\examples{
ma <- cbind(1,1,c(1,1,7))
D <- t(ma \%*\% c(0,1,-1))
## IGNORE_RDIFF_BEGIN
## note that results may vary according to BLAS
isKerAinKerB(D,ma)
isKerAinKerB(ma,D)
## IGNORE_RDIFF_END
}
\keyword{algebra}
\keyword{array}
\concept{pseudo inverse}
\concept{projector}
\concept{ker}
