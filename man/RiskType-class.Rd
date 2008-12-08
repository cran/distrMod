\name{RiskType-class}
\docType{class}
\alias{RiskType-class}
\alias{show,RiskType-method}
\alias{type,RiskType-method}

\title{Risk}
\description{Class of risks; e.g., estimator risks.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      type of risk. }
  }
}
\section{Methods}{
  \describe{
    \item{type}{\code{signature(object = "RiskType")}: 
      accessor function for slot \code{type}. }
    \item{show}{\code{signature(object = "RiskType")}}
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{}
%\examples{}
\concept{risk}
\keyword{classes}
