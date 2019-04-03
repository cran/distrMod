\name{L2ParamFamily-class}
\docType{class}
\alias{plot}
\alias{plot-methods}
\alias{L2ParamFamily-class}
\alias{FisherInfo}
\alias{FisherInfo,L2ParamFamily,missing-method}
\alias{FisherInfo,L2ParamFamily,ParamFamParameter-method}
\alias{L2deriv}
\alias{L2deriv,L2ParamFamily,missing-method}
\alias{L2deriv,L2ParamFamily,ParamFamParameter-method}
\alias{L2derivSymm}
\alias{L2derivSymm,L2ParamFamily-method}
\alias{L2derivDistr}
\alias{L2derivDistr,L2ParamFamily-method}
\alias{L2derivDistrSymm}
\alias{L2derivDistrSymm,L2ParamFamily-method}
\alias{checkL2deriv,L2ParamFamily-method}
\alias{E,L2ParamFamily,EuclRandVariable,missing-method}
\alias{E,L2ParamFamily,EuclRandMatrix,missing-method}
\alias{E,L2ParamFamily,EuclRandVarList,missing-method}
%\alias{modifyModel}
%\alias{modifyModel-methods}
%\alias{modifyModel,L2ParamFamily,ParamFamParameter-method}
\alias{plot,L2ParamFamily,missing-method}

\title{L2 differentiable parametric family}
\description{Class of L2 differentiable parametric families.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("L2ParamFamily", ...)}.
  More frequently they are created via the generating function 
  \code{L2ParamFamily}. 
}
\section{Slots}{
   \describe{
   \item{\code{name}}{[inherited from class \code{"ProbFamily"}] 
      object of class \code{"character"}: 
      name of the family. }
    \item{\code{distribution}}{[inherited from class \code{"ProbFamily"}] 
      object of class \code{"Distribution"}:
      member of the family. }
    \item{\code{distrSymm}}{[inherited from class \code{"ProbFamily"}] 
      object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
    \item{\code{param}}{[inherited from class \code{"ParamFamily"}]
      object of class \code{"ParamFamParameter"}:
      parameter of the family. }
    \item{\code{fam.call}}{[inherited from class \code{"ParamFamily"}]
      object of class \code{"call"}:
      call by which parametric family was produced.}
    \item{\code{makeOKPar}}{[inherited from class \code{"ParamFamily"}]
           object of class \code{"function"}:
           has argument \code{param} --- the (total) parameter, 
           returns valid parameter; used if \code{optim} resp. \code{optimize}---
           try to use ``illegal'' parameter values; then \code{makeOKPar} makes
           a valid parameter value out of the illegal one.}
    \item{\code{startPar}}{[inherited from class \code{"ParamFamily"}]
           object of class \code{"function"}:
           has argument \code{x} --- the data, 
           returns starting parameter for \code{optim} resp. \code{optimize}---
           a starting estimator in case parameter is multivariate
           or a search interval in case parameter is univariate.}
    \item{\code{modifyParam}}{[inherited from class \code{"ParamFamily"}]
      object of class \code{"function"}:
      mapping from the parameter space (represented by \code{"param"}) 
      to the distribution space (represented by \code{"distribution"}). }
    \item{\code{props}}{[inherited from class \code{"ProbFamily"}]
      object of class \code{"character"}:
      properties of the family. }
    \item{\code{L2deriv}}{
      object of class \code{"EuclRandVariable"}:
      L2 derivative of the family. Its map slot must contain a list of functions.
      Each function in this list must have just one argument \code{x},
      which is vectorized, (i.e., callable for a vector-valued \code{x}),
      and has a one-dimensional, numeric return value.}
    \item{\code{L2deriv.fct}}{
      object of class \code{"function"}: mapping from 
      the parameter space (argument \code{param} of class 
      \code{"ParamFamParameter"}) to a mapping from observation \code{x} to the  
      value of the L2derivative; \code{L2deriv.fct} is then used from observation
      \code{x} to value of the L2derivative; \code{L2deriv.fct} is used by 
      \code{modifyModel} to  move the L2deriv according to a change in the 
      parameter.
    More specifically, let us call the parts \code{main} and \code{nuisance}
    of the parameter the \emph{unknown} parameter. If this unknown parameter is
    one-dimensional, the return value of \code{L2deriv.fct} must be a function
    in argument \code{x}, which is vectorized, (i.e.,
    callable for a vector-valued \code{x}), and has a one-dimensional, numeric
    return value. In case the dimension of the unknown parameter is larger
    than one, the return value must be a list of functions, each of which
    satisfies the conditions formulated for the case of a one-dimensional
    parameter of interest. The order of the components of this list is
    the same as the order of the parameter coordinates in \code{main}, followed
    by the ones in \code{nuisance}.}
  \item{L2derivSymm}{ object of class \code{"FunSymmList"}:
    symmetry of the maps contained in \code{L2deriv}; a list
    of symmetry properties of the same length as the return value of
    \code{ L2deriv.fct }.}
  \item{L2derivDistr}{object of class \code{"OptionalDistrListOrCall"}
     (i.e., \code{NULL} or an object of class \code{"DistrList"} or
     the respective call to generate the latter object): if non-null
     and non-call, a list which includes the distribution of \code{L2deriv};
     the length of this list of univariate distributions must be of the same
     length as the return value of \code{ L2deriv.fct }.}
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}:
    symmetry of the distributions contained in \code{L2derivDistr};
    the length of this list of symmetry properties must be
    of the same length as the  return value of \code{ L2deriv.fct }. }
    \item{\code{FisherInfo.fct}}{
      object of class \code{"function"}: 
      mapping from the parameter space (argument  \code{param} of class 
      \code{"ParamFamParameter"}) to the set of positive
      semidefinite matrices; \code{FisherInfo.fct} is used by \code{modifyModel} to 
      move the Fisher information according to a change in the parameter } 
    \item{\code{FisherInfo}}{
      object of class \code{"PosDefSymmMatrix"}:
      Fisher information of the family. }
    \item{\code{.withEvalL2derivDistr}}{
      logical  of length one: if \code{TRUE} slot \code{L2derivDistr} gets
      evaluated, otherwise it is only kept as call. }
}}
\section{Extends}{
Class \code{"ParamFamily"}, directly.\cr
Class \code{"ProbFamily"}, by class \code{"ParamFamily"}.
}
\section{Methods}{
  \describe{
    \item{L2deriv}{\code{signature(object = "L2ParamFamily")}: 
      accessor function for \code{L2deriv}. }

    \item{L2deriv}{\code{signature(object = "L2ParamFamily", 
                                   param = "ParamFamParameter")}: 
      returns the L2derivative at \code{param}, i.e.
      evaluates slot function \code{L2deriv.fct} at \code{param}. }

    \item{L2derivSymm}{\code{signature(object = "L2ParamFamily")}: 
      accessor function for \code{L2derivSymm}. }

    \item{L2derivDistr}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{L2derivDistr}. }

    \item{L2derivDistrSymm}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{L2derivDistrSymm}. }

    \item{FisherInfo}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{FisherInfo}. }

    \item{FisherInfo}{\code{signature(object = "L2ParamFamily", 
                                   param = "ParamFamParameter")}: 
      returns the Fisher Information at \code{param}, i.e.
      evaluates slot function \code{FisherInfo.fct} at \code{param}. }

    \item{checkL2deriv}{\code{signature(object = "L2ParamFamily")}: 
      check centering of \code{L2deriv} and compute precision 
      of Fisher information. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandVariable", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandMatrix", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandVarList", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{plot}{\code{signature(x = "L2ParamFamily")}: 
      plot of \code{distribution} and \code{L2deriv}. More precisely,
      this method has arguments
      \code{
      plot(x, withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], ...,
             mfColRow = TRUE, to.draw.arg = NULL, withSubst = TRUE)
      } where \describe{
        \item{x}{object of class \code{"L2ParamFamily"}}
        \item{withSweave}{logical: if \code{TRUE} (for working with \command{Sweave}) 
                          no extra device is opened and height/width are not set}
        \item{main}{logical: is a main title to be used? or \cr
                    just as argument \code{main} in \code{\link[graphics]{plot.default}}.}
        \item{inner}{logical: do panels have their own titles? or \cr
               character vector of / cast to length 'number of plotted 
               panels' with the corresponding panel titles. For further 
               information, see also 
               \code{\link[distr:plot-methods]{plot}} and the 
               description of argument \code{main} in \code{\link[graphics]{plot.default}}.} 
        \item{sub}{logical: is a sub-title to be used? or \cr
                   just as argument \code{sub} in \code{\link[graphics]{plot.default}}.}
        \item{tmar}{top margin -- useful for non-standard main title sizes}
        \item{bmar}{bottom margin -- useful for non-standard sub title sizes}
        \item{cex.inner}{magnification to be used for inner titles relative
          to the current setting of \code{cex}; as in 
          \code{\link[graphics]{par}}; can be a vector of length 2; in this
          case the first component is for the distribution panels, the 
          second for the L2-derivative-panels.}
       \item{col.inner}{character or integer code; color for the inner title}              
       \item{mfColRow}{shall default partition in panels be used --- 
                       defaults to \code{TRUE}}
       \item{to.draw.arg}{Either \code{NULL} (default; 
                          everything is plotted) or a vector of either integers 
                         (the indices of the subplots to be drawn) or characters 
                         --- the names of the subplots to be drawn: these
                         names are to be chosen among 
                         \code{c("d","p","q", dimnms)} where \code{dimnms} is
                         either the row names of the trafo matrix 
                         \code{rownames(trafo(x@param))} or if the last expression
                         is \code{NULL} a vector \code{"dim<dimnr>"}, 
                         \code{dimnr} running through the number of rows of the 
                         trafo matrix.
                         }
       \item{withSubst}{logical; if \code{TRUE} (default) pattern substitution for
          titles and lables is used; otherwise no substitution is used. }

       \item{\dots}{addtional arguments for \code{plot} --- see 
               \code{\link[graphics]{plot}}, 
               \code{\link[graphics]{plot.default}}, 
               \code{\link[stats]{plot.stepfun}} }

      } 
      If \code{\dots} contains argument \code{ylim}, this may either be
      as in \code{plot.default} (i.e. a vector of length 2) or a vector of 
      length 4, where the first two elements are the values for \code{ylim} 
      in panels "d.c" and "d.d", and the last
      two elements are the values for \code{ylim} resp. \code{xlim} in panels 
      "p", "p.c", "p.d" and "q", "q.c", "q.d".
      In all title and axis label arguments, if \code{withSubst} is \code{TRUE},
      the following patterns are substituted:
      \describe{
         \item{\code{"\%C"}}{class of argument \code{x}}
         \item{\code{"\%A"}}{deparsed argument \code{x}}
         \item{\code{"\%D"}}{time/date-string when the plot was generated}
      }

      In addition, argument \code{\dots} may contain arguments \code{panel.first},
      \code{panel.last}, i.e., hook expressions to be evaluated at the very beginning
      and at the very end of each panel (within the then valid coordinates).
      To be able to use these hooks for each panel individually, they may also be
     lists of expressions (of the same length as the number of panels and
     run through in the same order as the panels).

     The return value of the plot methods is an S3 object of class
     \code{c("plotInfo","DiagnInfo")}, i.e., a list
     containing the information needed to produce the
     respective plot, which at a later stage could be used by different
     graphic engines (like, e.g. \code{ggplot}) to produce the plot
     in a different framework. A more detailed description will follow in
     a subsequent version.
}

    \item{modifyModel}{\code{signature(model = "L2ParamFamily", param = "ParamFamParameter")}:
      moves the L2-parametric Family \code{model} to parameter \code{param} }
  }
}

\details{
In the \code{E}-methods, diagnostics on the involved integrations are available
if argument   \code{diagnostic} is \code{TRUE}. Then there is attribute
\code{diagnostic} attached to the return value, which may be inspected
  and accessed through \code{\link[distrEx:distrExIntegrate]{showDiagnostic}} and
  \code{\link[distrEx:distrExIntegrate]{getDiagnostic}}.

}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily}}, \code{\link{ParamFamily-class}}}
\examples{
F1 <- new("L2ParamFamily")
plot(F1)

## selection of subpanels for plotting
F2 <- L2LocationScaleFamily()
layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE))
plot(F2,mfColRow = FALSE,
     to.draw.arg=c("p","q","loc"))
plot(F2,mfColRow = FALSE, inner=list("empirical cdf","pseudo-inverse",
     "L2-deriv, loc.part"), to.draw.arg=c("p","q","loc"))
}
\concept{parametric family}
\keyword{classes}
\keyword{models}
