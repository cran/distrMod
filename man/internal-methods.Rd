\name{.checkEstClassForParamFamily-methods}
\docType{methods}
\alias{.checkEstClassForParamFamily}
\alias{.checkEstClassForParamFamily-methods}
\alias{.checkEstClassForParamFamily,ANY-method}
\alias{.checkEstClassForParamFamily,ANY,ANY-method}
\title{ Methods for Function .checkEstClassForParamFamily in Package `distrMod' }

\description{.checkEstClassForParamFamily-methods}

\usage{
.checkEstClassForParamFamily(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{ANY,ANY}(PFam, estimator)
}

\arguments{
  \item{PFam}{a parametric family. }
  \item{estimator}{an estimator. }
  }

\details{
The respective methods can be used to cast an estimator to
a model-specific subclass with particular methods.
}

\value{
The (default) \code{ANY,ANY}-method returns the estimator unchanged.
}
\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  }

\keyword{methods}

