\name{NbinomFamily}
\alias{NbinomFamily}
\alias{NbinomwithSizeFamily}
\alias{NbinomMeanSizeFamily}

\title{Generating function for Nbinomial families}
\description{
  Generates an object of class \code{"L2ParamFamily"} which
  represents a Nbinomial family where the probability of
  success is the parameter of interest.
}
\usage{
NbinomFamily(size = 1, prob = 0.5, trafo)
NbinomwithSizeFamily(size = 1, prob = 0.5, trafo, withL2derivDistr = TRUE)
NbinomMeanSizeFamily(size = 1, mean = 0.5, trafo, withL2derivDistr = TRUE )
}
\arguments{
  \item{size}{ number of trials }
  \item{prob}{ probability of success }
  \item{mean}{ alternative parameter for negative binomial parameter}
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
  \item{withL2derivDistr}{logical: shall the distribution of the L2 derivative
      be computed? Defaults to \code{TRUE}; setting it to \code{FALSE}
      speeds up computations.}
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
  \code{NbinomFamily} assumes \code{size} to be known; while
  for \code{NbinomwithSizeFamily} it is a second (unknown) parameter;
  for \code{NbinomMeanSizeFamily} is like \code{NbinomwithSizeFamily}
  but uses the \code{size,mean} parametrization instead of the
  \code{size,prob} one.  
}
\value{Object of class \code{"L2ParamFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the 
  Asymptotic Theory of Robustness}. Bayreuth: Dissertation.

  Kohl, M. and Ruckdeschel, P. (2010). 
  R Package distrMod: S4 Classes and Methods for
  Probability Models. To appear in Journal of Statistical Software.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link[distr]{Nbinom-class}}}
\examples{
(N1 <- NbinomFamily(size = 25, prob = 0.25))
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)
(N1.w <- NbinomwithSizeFamily(size = 25, prob = 0.25))
plot(N1.w)
FisherInfo(N1.w)
checkL2deriv(N1.w)
(N2.w <- NbinomMeanSizeFamily(size = 25, mean = 75))
plot(N2.w)
FisherInfo(N2.w)
checkL2deriv(N2.w)

}
\concept{Negative Binomial model}
\keyword{models}
