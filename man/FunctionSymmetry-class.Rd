\name{FunctionSymmetry-class}
\docType{class}
\alias{FunctionSymmetry-class}
%\alias{coerce,FunctionSymmetry,FunSymmList-method}

\title{Class of Symmetries for Functions}
\description{Class of symmetries for functions.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      discribes type of symmetry. }
    \item{\code{SymmCenter}:}{Object of class \code{"OptionalNumeric"}:
      center of symmetry. }
  }
}
\section{Extends}{
Class \code{"Symmetry"}, directly.
}
%\section{Methods}{
%  \describe{
%    \item{coerce}{\code{signature(from = "FunctionSymmetry", to = "FunSymmList")}: 
%    create a \code{"FunSymmList"} object from a \code{"FunctionSymmetry"} object.  }
%  }
%}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{Symmetry-class}}, \code{\link{OptionalNumeric-class}}}
%\examples{}
\concept{symmetry}
\keyword{classes}
