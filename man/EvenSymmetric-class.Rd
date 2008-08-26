\name{EvenSymmetric-class}
\docType{class}
\alias{EvenSymmetric-class}

\title{Class for Even Functions}
\description{Class for even functions.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("EvenSymmetric")}.
  More frequently they are created via the generating function 
  \code{EvenSymmetric}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      contains \dQuote{even function} }
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
\seealso{\code{\link{EvenSymmetric}}, \code{\link{FunctionSymmetry-class}}}
\examples{
new("EvenSymmetric")
}
\concept{even function}
\concept{symmetry}
\keyword{classes}
