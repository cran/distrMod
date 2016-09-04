\name{L2LocationFamily}
\alias{L2LocationFamily}

\title{Generating function for L2LocationFamily-class}
\description{
  Generates an object of class \code{"L2LocationFamily"}.
}
\usage{
L2LocationFamily(loc = 0, name, centraldistribution = Norm(),
                 locname = "loc", modParam, LogDeriv,  
                 L2derivDistr.0, FisherInfo.0, distrSymm, L2derivSymm, 
                 L2derivDistrSymm, trafo, .returnClsName = NULL)
}
\arguments{
  \item{loc}{numeric: location parameter of the model. }
  \item{name}{character: name of the parametric family. }
  \item{centraldistribution}{object of class \code{"AbscontDistribution"};
    we assume from the beginning, that centraldistribution
    is symmetric about its median.}
  \item{modParam}{ optional function: mapping from the parameter space 
    (represented by \code{"param"}) to the distribution space 
    (represented by \code{"distribution"}). }
  \item{locname}{a character vector of length 1 containing the name
  of the location parameter}
  \item{LogDeriv}{function with argument \code{x}: the negative logarithmic 
    derivative of the density of the central distribution; if missing, it is
    determined numerically using numeric differentiation. }
  \item{L2derivDistr.0}{object of class \code{"UnivariateDistribution"}:
    distribution of the L2derivative at the central distribution }
  \item{FisherInfo.0}{object of class \code{"PosSemDefSymmMatrix"}:
    Fisher information of the model at the "standard" parameter value}
  \item{distrSymm}{ object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
  \item{L2derivSymm}{ object of class \code{"FunSymmList"}: 
    symmetry of the maps contained in \code{L2deriv} }
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}: 
    symmetry of the distributions contained in \code{L2derivDistr} }
  \item{trafo}{ matrix or  function in \code{param}: transformation of the parameter }
  \item{.returnClsName}{ the class name of the return value; by default this
    argument is \code{NULL} whereupon the return class will be 
    \code{L2LocationScaleFamily}; but, internally, this generating function is also 
    used to produce objects of class Classes \code{NormLocationFamily} and 
    \code{GumbelLocationFamily} (the latter in package \pkg{RobExtremes}.}
}
\details{
  If \code{name} is missing, the default 
  \dQuote{L2 location family} is used.
  The function \code{modParam} is optional. If it is missing, it is 
  constructed from \code{centraldistribution} using the location structure
  of the model.
  Slot \code{param} is filled accordingly with the argument 
  \code{trafo} passed to \code{L2LocationFamily}.
  In case \code{L2derivDistr.0} is missing, \code{L2derivDistr} is computed 
  via \code{imageDistr}, else \code{L2derivDistr} is assigned 
  \code{L2derivDistr.0}, coerced to \code{"UnivariateDistributionList"}.
  In case \code{FisherInfo.0} is missing, Fisher information is computed 
  from \code{L2deriv} using \code{E}.
  If \code{distrSymm} is missing, it is set to symmetry about \code{loc}.
  If \code{L2derivSymm} is missing, it is set to no symmetry, and
  if \code{L2derivDistrSymm} is missing, it is set to no symmetry, too. 
}
\value{Object of class \code{"L2LocationFamily"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{L2LocationFamily-class}}}
\examples{
F1 <- L2LocationFamily()
plot(F1)
}
\concept{parametric family}
\keyword{models}
