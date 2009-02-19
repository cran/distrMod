\name{existsPIC-methods}
\docType{methods}
\alias{existsPIC-methods}
\alias{existsPIC}
\alias{existsPIC,L2ParamFamily-method}
\title{ Methods for Function existsPIC in Package `distrMod'}
\description{existsPIC-methods to check whether in a given L2 differentiable model at parameter value theta there exist
(partial) influence curves to Trafo \eqn{D_\theta}{D_theta}.}

\usage{existsPIC(object, ...)
\S4method{existsPIC}{L2ParamFamily}(object, warning = TRUE, tol = .Machine$double.eps)
}
\arguments{
  \item{object}{L2ParamFamily}
  \item{...}{further arguments used by specific methods.}
  \item{warning}{logical: should a warning be issued if there exist no (partial) influence curves?}
  \item{tol}{ the tolerance the linear algebraic operations.  Default is \code{.Machine$double.eps}.}
}
\details{To check the existence of (partial) influence curves and, simultaneously, for bounded
(partial) influence curves, by Lemma 1.1.3 in Kohl(2005) [resp. the fact that \eqn{{\rm ker} I={\rm ker} J}{ker I = ker J}
for \eqn{J= {\rm E} (\Lambda',1)' (\Lambda',1) w}{J= E (Lambda',1)' (Lambda',1) w} and
\eqn{w={\rm min}(1, b/|(\Lambda',1)|}{w=min(1, b/|(Lambda',1)|)}], it suffices to check that
\eqn{{\rm ker }I}{ker I} is a subset of \eqn{{\rm ker }D_\theta}{ker D_theta}. This is done by
a call to \code{isKerAinKerB}.}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
\seealso{
\code{\link{isKerAinKerB}}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.
}
\concept{partial influence curve}
\concept{bounded influence curve}
\concept{existence of influence curves}
\concept{influence curve}
\keyword{robust}
