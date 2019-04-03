\name{L2ParamFamily}
\alias{L2ParamFamily}

\title{Generating function for L2ParamFamily-class}
\description{
  Generates an object of class \code{"L2ParamFamily"}.
}
\usage{
L2ParamFamily(name, distribution = Norm(), distrSymm, 
              main = main(param), nuisance = nuisance(param),
              fixed = fixed(param), trafo = trafo(param),
              param = ParamFamParameter(name = paste("Parameter of", name),  
                          main = main, nuisance = nuisance, 
                          fixed = fixed, trafo = trafo),
              props = character(0),
              startPar = NULL, makeOKPar = NULL,
              modifyParam = function(theta){ Norm(mean=theta) },
              L2deriv.fct = function(param) {force(theta <- param@main)
                           return(function(x) {x-theta})},
              L2derivSymm, L2derivDistr, L2derivDistrSymm, 
              FisherInfo.fct, FisherInfo = FisherInfo.fct(param),
              .returnClsName = NULL, .withMDE = TRUE)
}
\arguments{
  \item{name}{ character string: name of the family }
  \item{distribution}{ object of class \code{"Distribution"}:
    member of the family }
  \item{distrSymm}{ object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
  \item{main}{ numeric vector: main parameter }
  \item{nuisance}{ numeric vector: nuisance parameter }
  \item{fixed}{ numeric vector: fixed part of the parameter }
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
  \item{param}{ object of class \code{"ParamFamParameter"}:
    parameter of the family }
  \item{startPar}{\code{startPar} is a function in the observations \code{x}
    returning initial information for \code{MCEstimator} used 
    by \code{optimize} resp. \code{optim}; i.e; if (total) parameter is of 
    length 1, \code{startPar} returns a search interval, else it returns an initial 
    parameter value.}
  \item{makeOKPar}{\code{makeOKPar} is a function in the (total)
    parameter \code{param}; used if \code{optim} resp. \code{optimize}---
           try to use ``illegal'' parameter values; then \code{makeOKPar} makes
           a valid parameter value out of the illegal one; if \code{NULL}
    slot \code{makeOKPar} of \code{ParamFamily} is used to produce it.}
  \item{modifyParam}{ function: mapping from the parameter space 
  (represented by \code{"param"}) to the distribution space 
  (represented by \code{"distribution"}). }
  \item{props}{ character vector: properties of the family }
  \item{L2deriv.fct}{ function: mapping from the parameter space (argument
    \code{param} of class \code{"ParamFamParameter"}) to a mapping from 
    observation \code{x} to the value of the L2derivative;  
    \code{L2deriv.fct} is used by \code{modifyModel} to 
    move the L2deriv according to a change in the parameter,
    and to fill slot \code{L2deriv}.
    More specifically, let us call the parts \code{main} and \code{nuisance}
    of the parameter the \emph{unknown} parameter. If this unknown parameter is
    one-dimensional, the return value of \code{L2deriv.fct} must be a function
    in argument \code{x}, which is vectorized, (i.e.,
    callable for a vector-valued \code{x}), and has a one-dimensional, numeric
    return value. In case the dimension of the unknown parameter is larger
    than one, the return value must be a list of functions, each of which
    satisfies the conditions formulated for the case of a one-dimensional
    parameter of interest. The order of the components of this list is
    the same as the order of the parameter coordinates in \code{main}, followed
    by the ones in \code{nuisance}.}
  \item{L2derivSymm}{ object of class \code{"FunSymmList"}: 
    symmetry of the maps contained in \code{L2deriv}; a list
    of symmetry properties of the same length as the return value of
    \code{ L2deriv.fct }.}
  \item{L2derivDistr}{ object of class \code{"UnivarDistrList"}: 
    distribution of \code{L2deriv}; the length of this list
    of univariate distributions must be of the same length as the
    return value of \code{ L2deriv.fct }.}
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}: 
    symmetry of the distributions contained in \code{L2derivDistr};
    the length of this list of symmetry properties must be
    of the same length as the  return value of \code{ L2deriv.fct }. }
  \item{FisherInfo.fct}{function: mapping from the parameter space (argument
    \code{param} of class \code{"ParamFamParameter"}) to the set of positive
    semidefinite matrices; \code{FisherInfo.fct} is used by \code{modifyModel} to 
    move the Fisher information according to a change in the parameter } 
  \item{FisherInfo}{ object of class \code{"PosSemDefSymmMatrix"}:
    Fisher information of the family }
   \item{.returnClsName}{ the class name of the return value; by default this
   argument is \code{NULL} whereupon the return class will be 
   \code{L2ParamFamily}; but, internally, this generating function is also 
   used to e.g. produce objects of class \code{BinomialFamily}, \code{PoisFamily}
   \code{GammaFamily}, \code{BetaFamily}.}
    \item{.withMDE}{ logical of length 1: Tells R how to use the function from
    slot \code{startPar} in case of a \code{kStepEstimator}---use it as is or
    to compute the starting point for a minimum distance estimator which in
    turn then serves as starting point for \code{roptest} / \code{robest}
    (from package \pkg{ROptEst}). If \code{TRUE} (default) the latter
    alternative is used. Ignored if \pkg{ROptEst} is not used.}
}
\details{
  If \code{name} is missing, the default 
  \dQuote{L2 differentiable parametric family of probability measures} 
  is used. In case \code{distrSymm} is missing it is set to
  \code{NoSymmetry()}.
  If \code{param} is missing, the parameter is created via
  \code{main}, \code{nuisance} and \code{trafo} as described
  in \code{\link{ParamFamParameter}}. In case \code{L2derivSymm} is
  missing, it is filled with an object of class \code{FunSymmList} 
  with entries \code{NonSymmetric()}. In case \code{L2derivDistr} is missing, 
  it is computed via \code{imageDistr}. If \code{L2derivDistrSymm} is missing, 
  it is set to an object of class \code{DistrSymmList} with entries 
  \code{NoSymmetry()}. In case \code{FisherInfo} is missing, it is computed 
  from \code{L2deriv} using \code{E}.
}
\value{Object of class \code{"L2ParamFamily"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}}
\examples{
F1 <- L2ParamFamily()
plot(F1)
}
\concept{parametric family}
\keyword{models}
