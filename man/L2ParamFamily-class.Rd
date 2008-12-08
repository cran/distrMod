\name{L2ParamFamily-class}
\docType{class}
\alias{L2ParamFamily-class}
\alias{FisherInfo}
\alias{FisherInfo,L2ParamFamily,missing-method}
\alias{FisherInfo,L2ParamFamily,ParamFamParameter-method}
\alias{L2deriv}
\alias{L2deriv,L2ParamFamily,missing-method}
\alias{L2deriv,L2ParamFamily,ParamFamParameter-method}
\alias{L2derivSymm}
\alias{L2derivSymm,L2ParamFamily-method}
\alias{L2derivDistr}
\alias{L2derivDistr,L2ParamFamily-method}
\alias{L2derivDistrSymm}
\alias{L2derivDistrSymm,L2ParamFamily-method}
\alias{checkL2deriv,L2ParamFamily-method}
\alias{E,L2ParamFamily,EuclRandVariable,missing-method}
\alias{E,L2ParamFamily,EuclRandMatrix,missing-method}
\alias{E,L2ParamFamily,EuclRandVarList,missing-method}
%\alias{modifyModel}
%\alias{modifyModel-methods}
%\alias{modifyModel,L2ParamFamily,ParamFamParameter-method}
\alias{plot,L2ParamFamily,missing-method}

\title{L2 differentiable parametric family}
\description{Class of L2 differentiable parametric families.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("L2ParamFamily", ...)}.
  More frequently they are created via the generating function 
  \code{L2ParamFamily}. 
}
\section{Slots}{
   \describe{
   \item{\code{name}:}{[inherited from class \code{"ProbFamily"}] 
      object of class \code{"character"}: 
      name of the family. }
    \item{\code{distribution}:}{[inherited from class \code{"ProbFamily"}] 
      object of class \code{"Distribution"}:
      member of the family. }
    \item{\code{distrSymm}:}{[inherited from class \code{"ProbFamily"}] 
      object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
    \item{\code{param}:}{[inherited from class \code{"ParamFamily"}]
      object of class \code{"ParamFamParameter"}:
      parameter of the family. }
    \item{\code{fam.call}:}{[inherited from class \code{"ParamFamily"}]
      object of class \code{"call"}:
      call by which parametric family was produced.}
    \item{\code{makeOKPar}:}{[inherited from class \code{"ParamFamily"}]
           object of class \code{"function"}:
           has argument \code{param} --- the (total) parameter, 
           returns valid parameter; used if \code{optim} resp. \code{optimize}---
           try to use ``illegal'' parameter values; then \code{makeOKPar} makes
           a valid parameter value out of the illegal one.}
    \item{\code{startPar}:}{[inherited from class \code{"ParamFamily"}]
           object of class \code{"function"}:
           has argument \code{x} --- the data, 
           returns starting parameter for \code{optim} resp. \code{optimize}---
           a starting estimator in case parameter is multivariate
           or a search interval in case parameter is univariate.}
    \item{\code{modifyParam}:}{[inherited from class \code{"ParamFamily"}]
      object of class \code{"function"}:
      mapping from the parameter space (represented by \code{"param"}) 
      to the distribution space (represented by \code{"distribution"}). }
    \item{\code{props}:}{[inherited from class \code{"ProbFamily"}]
      object of class \code{"character"}:
      properties of the family. }
    \item{\code{L2deriv}:}{
      object of class \code{"EuclRandVariable"}:
      L2 derivative of the family. }
    \item{\code{L2deriv.fct}:}{
      object of class \code{"function"}: mapping from 
      the parameter space (argument \code{param} of class 
      \code{"ParamFamParameter"}) to a mapping from observation \code{x} to the  
      value of the L2derivative; \code{L2deriv.fct} is then used from observation
      \code{x} to value of the L2derivative; \code{L2deriv.fct} is used by 
      \code{modifyModel} to  move the L2deriv according to a change in the 
      parameter }
    \item{\code{L2derivSymm}:}{[
      object of class \code{"FunSymmList"}:
      symmetry of the maps included in \code{L2deriv}. }
    \item{\code{L2derivDistr}:}{
      object of class \code{"UnivarDistrList"}:
      list which includes the distribution of \code{L2deriv}. }
    \item{\code{L2derivDistrSymm}:}{
      object of class \code{"DistrSymmList"}:
      symmetry of the distributions included in \code{L2derivDistr}. }
    \item{\code{FisherInfo.fct}:}{
      object of class \code{"function"}: 
      mapping from the parameter space (argument  \code{param} of class 
      \code{"ParamFamParameter"}) to the set of positive
      semidefinite matrices; \code{FisherInfo.fct} is used by \code{modifyModel} to 
      move the Fisher information according to a change in the parameter } 
    \item{\code{FisherInfo}:}{
      object of class \code{"PosDefSymmMatrix"}:
      Fisher information of the family. }
}}
\section{Extends}{
Class \code{"ParamFamily"}, directly.\cr
Class \code{"ProbFamily"}, by class \code{"ParamFamily"}.
}
\section{Methods}{
  \describe{
    \item{L2deriv}{\code{signature(object = "L2ParamFamily")}: 
      accessor function for \code{L2deriv}. }

    \item{L2deriv}{\code{signature(object = "L2ParamFamily", 
                                   param = "ParamFamParameter")}: 
      returns the L2derivative at \code{param}, i.e.
      evaluates slot function \code{L2deriv.fct} at \code{param}. }

    \item{L2derivSymm}{\code{signature(object = "L2ParamFamily")}: 
      accessor function for \code{L2derivSymm}. }

    \item{L2derivDistr}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{L2derivDistr}. }

    \item{L2derivDistrSymm}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{L2derivDistrSymm}. }

    \item{FisherInfo}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{FisherInfo}. }

    \item{FisherInfo}{\code{signature(object = "L2ParamFamily", 
                                   param = "ParamFamParameter")}: 
      returns the Fisher Information at \code{param}, i.e.
      evaluates slot function \code{FisherInfo.fct} at \code{param}. }

    \item{checkL2deriv}{\code{signature(object = "L2ParamFamily")}: 
      check centering of \code{L2deriv} and compute precision 
      of Fisher information. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandVariable", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandMatrix", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandVarList", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{plot}{\code{signature(x = "L2ParamFamily")}: 
      plot of \code{distribution} and \code{L2deriv}. }

    \item{modifyModel}{\code{signature(model = "L2ParamFamily", param = "ParamFamParameter")}:
      moves the L2-parametric Family \code{model} to parameter \code{param} }
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily}}, \code{\link{ParamFamily-class}}}
\examples{
F1 <- new("L2ParamFamily")
plot(F1)
}
\concept{parametric family}
\keyword{classes}
\keyword{models}
