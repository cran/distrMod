\name{modifyModel-methods}
\docType{methods}
\alias{modifyModel-methods}
\alias{modifyModel}
\alias{modifyModel,ParamFamily,ParamFamParameter-method}
\alias{modifyModel,L2ParamFamily,ParamFamParameter-method}
\alias{modifyModel,L2LocationFamily,ParamFamParameter-method}
\alias{modifyModel,L2ScaleFamily,ParamFamParameter-method}
\alias{modifyModel,L2LocationScaleFamily,ParamFamParameter-method}
\alias{modifyModel,GammaFamily,ParamFamParameter-method}
\alias{modifyModel,ExpScaleFamily,ParamFamParameter-method}
\title{ Methods for function modifyModel in Package `distrMod' }
\description{
  Methods for function \code{modifyModel} in package \pkg{distrMod};
  \code{modifyModel} moves a model from one parameter value to
  another.
 }
\usage{
modifyModel(model, param,...)
\S4method{modifyModel}{ParamFamily,ParamFamParameter}(model,param, 
                       .withCall = TRUE, ...)
\S4method{modifyModel}{L2ParamFamily,ParamFamParameter}(model,param, 
                       .withCall = TRUE, .withL2derivDistr = TRUE, ...)
\S4method{modifyModel}{L2LocationFamily,ParamFamParameter}(model,param, ...)
\S4method{modifyModel}{L2ScaleFamily,ParamFamParameter}(model,param, ...)
\S4method{modifyModel}{L2LocationScaleFamily,ParamFamParameter}(model,
                       param, ...)
\S4method{modifyModel}{GammaFamily,ParamFamParameter}(model,param, ...)
\S4method{modifyModel}{ExpScaleFamily,ParamFamParameter}(model,param, ...)
}
\arguments{
  \item{model}{an object of class \code{ParamFamily}  --- the model to move.}
  \item{param}{an object of class \code{ParamFamParameter} --- the parameter to move to.}
  \item{.withCall}{logical: shall slot \code{fam.call} be updated?}
  \item{.withL2derivDistr}{logical: shall slot \code{L2derivDistr} be updated?}
  \item{\dots}{additional argument(s) for methods; not used so far}
}
\value{
  a corresponding instance of the model in argument \code{model} with moved
  parameters.
  }

\details{
  \code{modifyModel} is merely used internally for moving
  the model along modified parameter values during a model
  fit.
    
  It generally simply copies the original model and only
  modifies the affected slots, i.e. 
 \code{distribution}, the distribution of the observations,
 \code{param}, the parameter,
 \code{L2deriv}, the L2-derivative at the parameter,
 \code{L2FisherInfo}, the Fisher information at the parameter,
 the symmetry slots \code{distrSymm}, \code{L2derivSymm}, and 
 \code{L2derivDistrSymm},
 and, finally, \code{L2derivDistr} the (marginal) distribution(s)
 of the L2derivative. By default, also slot  \code{fam.call} is updated.

  In case \code{model} is of class \code{L2LocationFamily},
  \code{L2ScaleFamily}, or \code{L2LocationScaleFamily},
  symmetry slots are updated to be centered about the 
  median of the (central) distribution (assuming the latter is
  symmetric about the median); as an intermediate step, these methods
  call the general \code{modifyModel}-method for signature \code{L2ParamFamily};
  in this call, however, slot \code{fam.call} is not updated (this is
  the reason for argument \code{.withCall}); this is then
  done in the individual parts of the corresponding method.
  }

\keyword{models}
