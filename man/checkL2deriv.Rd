\name{checkL2deriv}
\alias{checkL2deriv}

\title{Generic function for checking L2-derivatives}
\description{
  Generic function for checking the L2-derivative of
  an L2-differentiable family of probability measures.
}
\usage{
checkL2deriv(L2Fam, ...)
}
\arguments{
  \item{L2Fam}{ L2-differentiable family of probability measures }
  \item{\dots}{ additional parameters }
}
\details{ 
  The precisions of the centering and the Fisher information
  are computed.
}
\value{The maximum deviation is returned.}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}}
\examples{
F1 <- new("L2ParamFamily")
checkL2deriv(F1)
}
\concept{parametric family}
\keyword{models}
