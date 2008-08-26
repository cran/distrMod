\name{FunSymmList-class}
\docType{class}
\alias{FunSymmList-class}

\title{List of Symmetries for a List of Functions}
\description{Create a list of symmetries for a list of functions}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("FunSymmList", ...)}.
  More frequently they are created via the generating function 
  \code{FunSymmList}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{ Object of class \code{"list"}. A list 
      of objects of class \code{"FunctionSymmetry"}. }
  }
}
\section{Extends}{
Class \code{"list"}, from data part.\cr
Class \code{"vector"}, by class \code{"list"}.
}
%\section{Methods}{
%  \describe{
%    \item{coerce}{\code{signature(from = "FunctionSymmetry", to = "FunSymmList")}: 
%      create a \newline \code{"FunSymmList"} object from a 
%      \code{"FunctionSymmetry"} object. }
%  }
%}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{FunctionSymmetry-class}}}
\examples{
new("FunSymmList", list(NonSymmetric(), EvenSymmetric(SymmCenter = 1), 
                        OddSymmetric(SymmCenter = 2)))
}
\concept{symmetry}
\keyword{classes}
