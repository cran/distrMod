\name{L2ScaleUnknownLocationFamily}
\alias{L2ScaleUnknownLocationFamily}

\title{Generating function for L2LocationScaleFamily-class in nuisance situation}
\description{
  Generates an object of class \code{"L2LocationScaleFamily"} 
  in the situation where scale is main, location nuisance parameter.
}
\usage{
L2ScaleUnknownLocationFamily(loc = 0, scale = 1, name, centraldistribution = Norm(),
                      locscalename = c("loc", "scale"), modParam, LogDeriv,  
                      L2derivDistr.0, FisherInfo.0, distrSymm, L2derivSymm, 
                      L2derivDistrSymm, trafo, .returnClsName = NULL)

}
\arguments{
  \item{loc}{numeric: location parameter of the model. }
  \item{scale}{positive number: scale of the model. }
  \item{name}{character: name of the parametric family. }
  \item{centraldistribution}{object of class \code{"AbscontDistribution"}:
    central distribution; we assume by default, that centraldistribution
    is symmetric about \eqn{0}}
  \item{modParam}{ optional function: mapping from the parameter space 
    (represented by \code{"param"}) to the distribution space 
    (represented by \code{"distribution"}). }
  \item{locscalename}{a character vector of length 2 containing the names
  of the location and scale parameter; either unnamed, then order must
  be \code{c(loc,scale)}, or named, then names must be \code{"loc"} and 
  \code{"scale"}}
  \item{LogDeriv}{function with argument \code{x}: the negative logarithmic 
    derivative of the density of the central distribution; if missing, it is
    determined numerically using numeric differentiation. }
  \item{L2derivDistr.0}{list of length 2 of objects of class \code{"UnivariateDistribution"}:
    (marginal) distributions of the coordinates of the L2derivative at the 
    central distribution }
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
    used to produce objects of class \code{NormalLocationScaleFamily}, 
    \code{CauchyLocationScaleFamily}.} 
}
\details{
  If \code{name} is missing, the default 
  \dQuote{L2 scale family with unknown location (as nuisance)} is used.  
  The function \code{modParam} is optional. If it is missing, it is 
  constructed from \code{centraldistribution} using the location and
  scale structure of the model.
  Slot \code{param} is filled accordingly with the argument 
  \code{trafo} passed to \code{L2ScaleUnknownLocationFamily}.
  In case \code{L2derivDistr.0} is missing, \code{L2derivDistr} is computed 
  via \code{imageDistr}, else \code{L2derivDistr} is assigned 
  \code{L2derivDistr.0}, coerced to \code{"UnivariateDistributionList"}.
  In case \code{FisherInfo.0} is missing, Fisher information is computed 
  from \code{L2deriv} using \code{E}.
  If \code{distrSymm} is missing, it is set to symmetry about \code{loc}.
  If \code{L2derivSymm} is missing, its location and scale components are set 
  to no symmetry, respectively.
  if \code{L2derivDistrSymm} is missing, its location and scale components are set 
  to no symmetry, respectively.
}
\value{Object of class \code{"L2LocationScaleFamily"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{L2LocationScaleFamily-class}}}
\examples{
F1 <- L2ScaleUnknownLocationFamily()
plot(F1)
}
\concept{parametric family}
\keyword{models}
