\name{ParamFamParameter-class}
\docType{class}
\alias{ParamFamParameter-class}
\alias{ParamWithScaleFamParameter-class}
\alias{ParamWithScaleAndShapeFamParameter-class}
\alias{ParamWithShapeFamParameter-class}
\alias{length,ParamFamParameter-method}
\alias{dimension,ParamFamParameter-method}
\alias{main}
\alias{main,ParamFamParameter-method}
\alias{main,ParamWithScaleAndShapeFamParameter-method}
\alias{main<-}
\alias{main<-,ParamFamParameter-method}
\alias{nuisance}
\alias{nuisance,ParamFamParameter-method}
\alias{nuisance,ParamWithScaleAndShapeFamParameter-method}
\alias{nuisance<-}
\alias{nuisance<-,ParamFamParameter-method}
\alias{fixed}
\alias{fixed,ParamFamParameter-method}
\alias{fixed,ParamWithScaleAndShapeFamParameter-method}
\alias{fixed<-}
\alias{fixed<-,ParamFamParameter-method}
\alias{withPosRestr}
\alias{withPosRestr,ParamWithShapeFamParameter-method}
\alias{withPosRestr<-}
\alias{withPosRestr<-,ParamWithShapeFamParameter-method}
\alias{show,ParamFamParameter-method}
\alias{show,ParamWithShapeFamParameter-method}
\alias{show,ParamWithScaleAndShapeFamParameter-method}

%\alias{trafo}
%\alias{trafo,ParamFamParameter-method}
%\alias{trafo<-}
%\alias{trafo<-,ParamFamParameter-method}

\title{Parameter of a parametric family of probability measures}
\description{Class of the parameter of parametric families 
  of probability measures.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ParamFamParameter", ...)}.
  More frequently they are created via the generating function 
  \code{ParamFamParameter}. 
}
\section{Slots}{
  \describe{
    \item{\code{main}}{Object of class \code{"numeric"}: main parameter. }
    \item{\code{nuisance}}{Object of class \code{"OptionalNumeric"}:
      optional nuisance parameter. }
    \item{\code{fixed}}{Object of class \code{"OptionalNumeric"}:
      optional fixed part of the parameter. }
    \item{\code{trafo}}{Object of class \code{"MatrixorFunction"}: 
      transformation of the parameter. }
    \item{\code{name}}{Object of class \code{"character"}:
      name of the parameter. }
    \item{\code{withPosRestr}}{(for \code{ParamWithShapeFamParameter} and
     \code{ParamWithScaleAndShapeFamParameter}): Object of class \code{"logical"}:
     Is shape restricted to be positive? }
  }
}
\section{Extends}{
Class \code{"Parameter"}, directly.\cr
Class \code{"OptionalParameter"}, by class \code{"Parameter"}.
}
\section{Methods}{
  \describe{
    \item{main}{\code{signature(object = "ParamFamParameter")}: 
      accessor function for slot \code{main}. }

    \item{main<-}{\code{signature(object = "ParamFamParameter")}: 
      replacement function for slot \code{main}. }

    \item{nuisance}{\code{signature(object = "ParamFamParameter")}: 
      accessor function for slot \code{nuisance}. }

    \item{nuisance<-}{\code{signature(object = "ParamFamParameter")}: 
      replacement function for slot \code{nuisance}. }

    \item{fixed}{\code{signature(object = "ParamFamParameter")}: 
      accessor function for slot \code{fixed}. }

    \item{fixed<-}{\code{signature(object = "ParamFamParameter")}: 
      replacement function for slot \code{fixed}. }

    \item{trafo}{\code{signature(object = "ParamFamParameter")}: 
      accessor function for slot \code{trafo}. }

    \item{trafo<-}{\code{signature(object = "ParamFamParameter")}: 
      replacement function for slot \code{trafo}. }

    \item{length}{\code{signature(x = "ParamFamParameter")}: 
      sum of the lengths of \code{main} and \code{nuisance}. }

    \item{dimension}{\code{signature(x = "ParamFamParameter")}: 
      length of \code{main}. }

    \item{withPosRestr}{\code{signature(object = "ParamWithShapeFamParameter")}:
      accessor function for slot \code{trafo}. }

    \item{withPosRestr<-}{\code{signature(object = "ParamWithShapeFamParameter")}:
      replacement function for slot \code{trafo}. }

    \item{show}{\code{signature(object = "ParamFamParameter")}}
    \item{show}{\code{signature(object = "ParamWithShapeFamParameter")}}
    \item{show}{\code{signature(object = "ParamWithScaleAndShapeFamParameter")}}
  }
}
\section{Details for methods 'show', 'print'}{
Detailedness of output by methods \code{show}, \code{print} is controlled
by the global option \code{show.details} to be set by 
\code{\link{distrModoptions}}.

As method \code{show} is used when inspecting an object by typing the object's 
name into the console, \code{show} comes without extra arguments and hence 
detailedness must be controlled  by global options. 

Method \code{print} may be called with a (partially matched) argument
\code{show.details}, and then the global option is temporarily set to this
value.

More specifically, when \code{show.detail} is matched to \code{"minimal"}
only class and name as well as main and nuisance part of the
parameter are shown. When \code{show.detail} is matched to \code{"medium"}, 
and if you estimate non-trivial (i.e. not the identity)
transformation of the parameter of the parametric family, you will in addition 
be shown the derivative matrix, if the transformation is given in form of
this matrix, while, if the transformation is in function form, you will
only be told this.
Finally, when \code{show.detail} is matched to \code{"maximal"}, 
and you have a non-trivial transformation in function form, you will
also be shown the code to this function.
            }
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[distr]{Parameter-class}}}
\examples{
new("ParamFamParameter")
}
\concept{parameter}
\keyword{classes}
