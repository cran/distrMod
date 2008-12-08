\name{InternalReturnClasses-class}
\docType{class}
\alias{InternalReturnClasses-class}
\alias{InternalReturnClasses-class}
\alias{BinomFamily-class}
\alias{PoisFamily-class}
\alias{NormLocationFamily-class}
\alias{GumbelLocationFamily-class}
\alias{NormScaleFamily-class}
\alias{ExpScaleFamily-class}
\alias{LnormScaleFamily-class}
\alias{GammaFamily-class}
\alias{BetaFamily-class}, 
\alias{NormLocationScaleFamily-class}
\alias{CauchyLocationScaleFamily-class}

\title{Internal return classes for generating functions}
\description{internal return classes for generating functions 'L2ParamFamily' and 
'L2LocationFamily' (and friends); used for particular method dispatch only}

\section{Described classes}{
In this file we describe classes 
\code{BinomFamily}, \code{PoisFamily}, \code{GammaFamily}
\code{BetaFamily} ``extending'' (no new slots!) class \code{L2ParamFamily},
classes \code{NormLocationFamily} and \code{GumbelLocationFamily}, 
``extending'' (no new slots!) class \code{"L2LocationFamily"}, classes
\code{NormScaleFamily}, \code{ExpScaleFamily}, and \code{LnormScaleFamily} 
``extending'' (no new slots!) class \code{"L2ScaleFamily"}, and classes
\code{CauchyLocationScaleFamily} and \code{NormLocationScaleFamily},
``extending'' (no new slots!) class \code{"L2LocationScaleFamily"}.
}
\section{Objects from these classes}{
  Objects are only generated internally by the mentioned generating functions. 
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
    \item{\code{L2deriv}:}{[inherited from class \code{"L2ParamFamily"}]
      object of class \code{"EuclRandVariable"}:
      L2 derivative of the family. }
    \item{\code{L2deriv.fct}:}{[inherited from class \code{"L2ParamFamily"}]
      object of class \code{"function"}: mapping from 
      the parameter space (argument \code{param} of class 
      \code{"ParamFamParameter"}) to a mapping from observation \code{x} to the  
      value of the L2derivative; \code{L2deriv.fct} is then used from observation
      \code{x} to value of the L2derivative; \code{L2deriv.fct} is used by 
      \code{modifyModel} to  move the L2deriv according to a change in the 
      parameter }
    \item{\code{L2derivSymm}:}{[inherited from class \code{"L2ParamFamily"}]
      object of class \code{"FunSymmList"}:
      symmetry of the maps included in \code{L2deriv}. }
    \item{\code{L2derivDistr}:}{[inherited from class \code{"L2ParamFamily"}]
      object of class \code{"UnivarDistrList"}:
      list which includes the distribution of \code{L2deriv}. }
    \item{\code{L2derivDistrSymm}:}{[inherited from class \code{"L2ParamFamily"}]
      object of class \code{"DistrSymmList"}:
      symmetry of the distributions included in \code{L2derivDistr}. }
    \item{\code{FisherInfo.fct}:}{[inherited from class \code{"L2ParamFamily"}]
      object of class \code{"function"}: 
      mapping from the parameter space (argument  \code{param} of class 
      \code{"ParamFamParameter"}) to the set of positive
      semidefinite matrices; \code{FisherInfo.fct} is used by \code{modifyModel} to 
      move the Fisher information according to a change in the parameter } 
    \item{\code{FisherInfo}:}{[inherited from class \code{"L2ParamFamily"}]
      object of class \code{"PosDefSymmMatrix"}:
      Fisher information of the family. }
    \item{\code{LogDeriv}:}{(only loc/scale classes)[inherited from class \code{"L2GroupParamFamily"}]
      object of class \code{"function"}: has argument \code{x}; 
      the negative logarithmic derivative of the density of the model distribution at 
      the "standard" parameter value. }
    \item{\code{locscalename}:}{(only loc/scale classes)[inherited from class \code{"L2LocationScaleUnion"}]
      object of class \code{"character"}: names of location and scale parameter }
}}



\section{Extends}{
Classes \code{BinomFamily}, \code{PoisFamily}, \code{GammaFamily}
\code{BetaFamily} ``extend'' (no new slots!):\cr
Class \code{"L2ParamFamily"},  directly.\cr
Class \code{"ParamFamily"}, by class \code{"L2ParamFamily"}.\cr
Class \code{"ProbFamily"}, by class \code{"ParamFamily"}.
{ }\cr
Classes \code{NormLocationFamily} and \code{GumbelLocationFamily}, 
``extend'' (no new slots!):\cr
Class \code{"L2LocationFamily"}, directly.\cr
Class \code{"L2LocationScaleUnion"}, by class \code{"L2LocationFamily"}.\cr
Class \code{"L2GroupParamFamily"}, by class \code{"L2LocationScaleUnion"}.\cr
Class \code{"L2ParamFamily"},  directly.\cr
Class \code{"ParamFamily"}, by class \code{"L2ParamFamily"}.\cr
Class \code{"ProbFamily"}, by class \code{"ParamFamily"}.
{ }\cr
\code{NormScaleFamily}, \code{ExpScaleFamily}, and \code{LnormScaleFamily} 
``extend'' (no new slots!):\cr
Class \code{"L2ScaleFamily"}, directly.\cr
Class \code{"L2LocationScaleUnion"}, by class \code{"L2ScaleFamily"}.\cr
Class \code{"L2GroupParamFamily"}, by class \code{"L2LocationScaleUnion"}.\cr
Class \code{"L2ParamFamily"},  directly.\cr
Class \code{"ParamFamily"}, by class \code{"L2ParamFamily"}.\cr
Class \code{"ProbFamily"}, by class \code{"ParamFamily"}.
{ }\cr
\code{CauchyLocationScaleFamily} and \code{NormLocationScaleFamily}
``extend'' (no new slots!):\cr
Class \code{"L2LocationScaleFamily"}, directly.\cr
Class \code{"L2LocationScaleUnion"}, by class \code{"L2LocationScaleFamily"}.\cr
Class \code{"L2GroupParamFamily"}, by class \code{"L2LocationScaleUnion"}.\cr
Class \code{"L2ParamFamily"},  directly.\cr
Class \code{"ParamFamily"}, by class \code{"L2ParamFamily"}.\cr
Class \code{"ProbFamily"}, by class \code{"ParamFamily"}.
}
\section{Methods}{
   not yet done...
}
\seealso{\code{\link[methods:BasicClasses]{numeric-class}},
\code{\link{L2ParamFamily-class}},
\code{\link{L2GroupParamFamily-class}},
\code{\link{L2LocationFamily-class}},
\code{\link{L2ScaleFamily-class}},
\code{\link{L2LocationScaleFamily-class}},
}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\concept{parametric family}
\keyword{classes}
\keyword{models}
\keyword{internal}
