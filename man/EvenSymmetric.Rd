\name{EvenSymmetric}
\alias{EvenSymmetric}

\title{Generating function for EvenSymmetric-class}
\description{Generates an object of class \code{"EvenSymmetric"}.}
\usage{EvenSymmetric(SymmCenter = 0)}
\arguments{
  \item{SymmCenter}{ numeric: center of symmetry }
}
%\details{}
\value{Object of class \code{"EvenSymmetric"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{EvenSymmetric-class}}, \code{\link{FunctionSymmetry-class}}}
\examples{
EvenSymmetric()

## The function is currently defined as
function(SymmCenter = 0){ 
    new("EvenSymmetric", SymmCenter = SymmCenter) 
}
}
\concept{even function}
\concept{symmetry}
\keyword{math}
