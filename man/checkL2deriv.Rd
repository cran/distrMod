\name{checkL2deriv}
\alias{checkL2deriv}
\alias{print.relMatrix}

\title{Generic function for checking L2-derivatives}
\description{
  Generic function for checking the L2-derivative of
  an L2-differentiable family of probability measures.
}
\usage{
checkL2deriv(L2Fam, ...)
\method{print}{relMatrix}(x,...)
}
\arguments{
  \item{L2Fam}{ L2-differentiable family of probability measures }
  \item{x}{ argument to be printed }
  \item{\dots}{ additional parameters (ignored/for compatibility with S3 generic 
                in case \code{print.relMatrix}) }
}
\details{ 
  The precisions of the centering and the Fisher information
  are computed.
}
\value{A list with items \code{maximum.deviation}, \code{cent}, \code{consist},
       and \code{condition} is invisibly returned, where \code{maximum.deviation}
       comprises the maximal absolute value of all entries in \code{cent} and \code{consist},
       \code{cent} shows the expectation of \code{L2deriv(L2Fam)}
       (which should be 0), \code{consist} shows the difference
       between the Fisher information and \code{cov(L2deriv(L2Fam))} (which should be 0), and
       \code{condition} is the condition number of the Fisher information.}

\note{The return value gives the non-rounded values (which will be machine dependent),
      whereas on argument \code{out==TRUE} (the default) we only issue the values up to
      5 digits which should be independent of the machine. For the output of relative
	  differences, we adjust accuracy to the size of the maximal (absolute) value of 
	  the Fisher information. In case of the consistency condition, at positions where
	  the denominator is 0, we print a "."; this is done through helper S3 method
	  \code{print.relMatrix}.}

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
