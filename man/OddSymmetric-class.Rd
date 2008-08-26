\name{OddSymmetric-class}
\docType{class}
\alias{OddSymmetric-class}

\title{Class for Odd Functions}
\description{Class for odd functions.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("OddSymmetric")}.
  More frequently they are created via the generating function 
  \code{OddSymmetric}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      contains \dQuote{odd function} }
    \item{\code{SymmCenter}:}{Object of class \code{"numeric"}:
      center of symmetry }
  }
}
\section{Extends}{
Class \code{"FunctionSymmetry"}, directly.\cr
Class \code{"Symmetry"}, by class \code{"FunctionSymmetry"}.
}
%\section{Methods}{}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{OddSymmetric}}, \code{\link{FunctionSymmetry-class}}}
\examples{
new("OddSymmetric")
}
\concept{odd function}
\concept{symmetry}
\keyword{classes}
