\name{InternalClassUnions}
\docType{class}
\alias{InternalClassUnions}
\alias{OptionalNumeric-class}
\alias{OptionalNumericOrMatrix-class}
\alias{MatrixorFunction-class}
\alias{ShowDetails-class}
\alias{L2LocationScaleUnion-class}
\alias{locscalename}
\alias{locscalename,L2LocationScaleUnion-method}
\alias{locscalename<-}
\alias{locscalename<-,L2LocationScaleUnion-method}

\title{Class unions in 'distrMod'}
\description{Class unions in package \pkg{distrMod} defined for internal 
 purposes; these are \code{OptionalNumeric}, \code{OptionalMatrix},
 \code{MatrixorFunction}, \code{ShowDetails}}

\section{Objects from the Class}{All of these classes are virtual: 
 No objects may be created from them.}
\details{
These classes are used internally to make available methods or to allow
slots of classes to be filled with varying types. In particular
\describe{
\item{\code{"OptionalNumeric"}}{may contain objects of class  \code{"numeric"}
      or \code{"NULL"}; it is used e.g. for slot \code{nuisance} of class
      \code{"ParamFamParameter"}, as it may or may not be present but if so
      it has to be numeric.}
\item{\code{"OptionalNumericOrMatrix"}}{may contain objects of class  
      \code{"OptionalNumeric"} or \code{"matrix"}; it is used e.g. for 
       slot \code{asvar} of class \code{"Estimate"}, as it may or may not be 
       present but if so it has to be a number (numeric) or a matrix.}
\item{\code{"MatrixorFunction"}}{may contain objects of class  
      \code{"OptionalFunction"} or \code{"matrix"}; it is used e.g. for 
       slot \code{trafo} of class \code{"ParamFamParameter"}, as it may or may not be 
       present and if it is present, it has to either be a function or a matrix,
       see \code{link{trafo-methods}}.}
\item{\code{"ShowDetails"}}{may contain objects of class  
      \code{"Estimate"}, \code{"MCEstimate"}, \code{"Confint"}, \code{"PosSemDefSymmMatrix"}
       \code{"ParamFamily"}, or \code{"ParamFamParameter"}; used to provide 
       sort of a ``\code{show} with extra arguments'', in form of a common 
       \code{print} method for these S4  classes, which essentially just 
       temporarily sets the global options according to the optional arguments
       \code{digits} and \code{show.details}, calls \code{show} and then re-sets
       the options to their global settings.}
\item{\code{"L2LocationScaleUnion"}}{is a proper, but virtual in-between class
      between class \code{"L2GroupParamFamily"} and \code{"L2LocationFamily"}, 
      \code{"L2ScaleFamily"}, and \code{"L2LocationScaleFamily"};
      in addition to class \code{"L2GroupParamFamily"} it has a slot
      \code{locscalename} (with corresponding accessor and replacement method)
      to capture the fact that location and scale 
      parameter may carry names other than \code{"loc"} resp. \code{"scale"};
      used to provide a common class for (parts of) methods \code{modifyModel}
      \code{locscalename} and \code{locscalename<-}.}
}
}

%\section{Methods}{}
%\references{}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[methods:BasicClasses]{numeric-class}},
\code{\link[methods:BasicClasses]{matrix-class}},
\code{\link[methods:BasicClasses]{function-class}},
\code{\link{MCEstimate-class}},
\code{\link{Estimate-class}},
\code{\link{Confint-class}},
\code{\link{PosSemDefSymmMatrix-class}},
\code{\link{L2GroupParamFamily-class}},
\code{\link{L2LocationFamily-class}},
\code{\link{L2ScaleFamily-class}},
\code{\link{L2LocationScaleFamily-class}},
\code{\link{ParamFamParameter-class}}
\code{\link{ParamFamily-class}}
}
\keyword{classes}
\keyword{internal}
