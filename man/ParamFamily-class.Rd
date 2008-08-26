\name{ParamFamily-class}
\docType{class}
\alias{ParamFamily-class}
\alias{main,ParamFamily-method}
\alias{nuisance,ParamFamily-method}
\alias{fixed,ParamFamily-method}
\alias{param,ParamFamily-method}
\alias{modifyParam}
\alias{modifyParam,ParamFamily-method}
\alias{fam.call}
\alias{fam.call,ParamFamily-method}
\alias{startPar}
\alias{startPar,ParamFamily-method}
\alias{makeOKPar}
\alias{makeOKPar,ParamFamily-method}
\alias{plot,ParamFamily,ANY-method}
\alias{show,ParamFamily-method}
%\alias{trafo,ParamFamily-method}

\title{Parametric family of probability measures.}
\description{Class of parametric families of probability measures.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ParamFamily", ...)}.
  More frequently they are created via the generating function 
  \code{ParamFamily}. 
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
    \item{\code{param}:}{
      object of class \code{"ParamFamParameter"}:
      parameter of the family. }
    \item{\code{fam.call}:}{
      object of class \code{"call"}:
      call by which parametric family was produced.}
    \item{\code{makeOKPar}:}{
           object of class \code{"function"}:
           has argument \code{param} --- the (total) parameter, 
           returns valid parameter; used if \code{optim} resp. \code{optimize}---
           try to use ``illegal'' parameter values; then \code{makeOKPar} makes
           a valid parameter value out of the illegal one.}
    \item{\code{startPar}:}{
           object of class \code{"function"}:
           has argument \code{x} --- the data, 
           returns starting parameter for \code{optim} resp. \code{optimize}---
           a starting estimator in case parameter is multivariate
           or a search interval in case parameter is univariate.}
    \item{\code{modifyParam}:}{
      object of class \code{"function"}:
      mapping from the parameter space (represented by \code{"param"}) 
      to the distribution space (represented by \code{"distribution"}). }
    \item{\code{props}:}{[inherited from class \code{"ProbFamily"}]
      object of class \code{"character"}:
      properties of the family. }
}}
\section{Extends}{
Class \code{"ProbFamily"}, directly.
}
\section{Methods}{
  \describe{
    \item{main}{\code{signature(object = "ParamFamily")}: 
      wrapped accessor function for slot \code{main} of
      slot \code{param}. }

    \item{nuisance}{\code{signature(object = "ParamFamily")}: 
      wrapped accessor function for slot \code{nuisance} 
      of slot \code{param}. }

    \item{fixed}{\code{signature(object = "ParamFamily")}: 
      wrapped accessor function for slot \code{fixed} 
      of slot \code{param}. }

    \item{trafo}{\code{signature(object = "ParamFamily", param = "missing")}: 
      wrapped accessor function for slot \code{trafo} 
      of slot \code{param}. }

    \item{param}{\code{signature(object = "ParamFamily")}: 
      accessor function for slot \code{param}. }

    \item{modifyParam}{\code{signature(object = "ParamFamily")}: 
      accessor function for slot \code{modifyParam}. }

    \item{fam.call}{\code{signature(object = "ParamFamily")}: 
      accessor function for slot \code{fam.call}. }

    \item{plot}{\code{signature(x = "ParamFamily")}: 
      plot of slot \code{distribution}. }

    \item{show}{\code{signature(object = "ParamFamily")}}

  }
}
\section{Details for methods \code{show}, \code{print}}{
Detailedness of output by methods \code{show}, \code{print} is controlled
by the global option \code{show.details} to be set by 
\code{\link{distrModoptions}}.

As method \code{show} is used when inspecting an object by typing the object's 
name into the console, \code{show} comes without extra arguments and hence 
detailedness must be controlled  by global options. 

Method \code{print} may be called with a (partially matched) argument
\code{show.details}, and then the global option is temporarily set to this
value.

For class \code{ParamFamily}, this becomes relevant for slot \code{param}.
For details therefore confer to \code{\link{ParamFamParameter-class}}.
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distr]{Distribution-class}}}
\examples{
F1 <- new("ParamFamily") # prototype
plot(F1)
}
\concept{parameteric family}
\keyword{classes}
\keyword{models}
