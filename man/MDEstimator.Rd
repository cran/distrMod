\name{MDEstimator}
\alias{MDEstimator}
\alias{CvMMDEstimator}
\alias{KolmogorovMDEstimator}
\alias{TotalVarMDEstimator}
\alias{HellingerMDEstimator}
\alias{CvMDist2}

\title{ Function to compute minimum distance estimates }
\description{
  The function \code{MDEstimator} provides a general way to compute 
  minimum distance estimates.
}
\usage{
MDEstimator(x, ParamFamily, distance = KolmogorovDist, dist.name, 
            paramDepDist = FALSE, startPar = NULL, Infos, trafo = NULL,
            penalty = 1e20, validity.check = TRUE, asvar.fct, na.rm = TRUE,
            ..., .withEvalAsVar = TRUE, nmsffx = "",
            .with.checkEstClassForParamFamily = TRUE)
CvMMDEstimator(x, ParamFamily, muDatOrMod = c("Mod","Dat", "Other"),
            mu = NULL, paramDepDist = FALSE, startPar = NULL, Infos,
            trafo = NULL, penalty = 1e20, validity.check = TRUE, 
            asvar.fct = .CvMMDCovariance, na.rm = TRUE, ...,
            .withEvalAsVar = TRUE, nmsffx = "",
            .with.checkEstClassForParamFamily = TRUE)
KolmogorovMDEstimator(x, ParamFamily, paramDepDist = FALSE, startPar = NULL, Infos, 
            trafo = NULL, penalty = 1e20, validity.check = TRUE, asvar.fct, 
            na.rm = TRUE, ..., .withEvalAsVar = TRUE, nmsffx = "",
            .with.checkEstClassForParamFamily = TRUE)
TotalVarMDEstimator(x, ParamFamily, paramDepDist = FALSE, startPar = NULL, Infos, 
            trafo = NULL, penalty = 1e20, validity.check = TRUE, asvar.fct, 
            na.rm = TRUE, ..., .withEvalAsVar = TRUE, nmsffx = "",
            .with.checkEstClassForParamFamily = TRUE)
HellingerMDEstimator(x, ParamFamily, paramDepDist = FALSE, startPar = NULL, Infos, 
            trafo = NULL, penalty = 1e20, validity.check = TRUE, asvar.fct, 
            na.rm = TRUE, ..., .withEvalAsVar = TRUE, nmsffx = "",
            .with.checkEstClassForParamFamily = TRUE)
CvMDist2(e1,e2,... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ (empirical) data }
  \item{ParamFamily}{ object of class \code{"ParamFamily"} }
  \item{distance}{ (generic) function: to compute distance beetween (emprical)
  data and objects of class \code{"Distribution"}. }
  \item{dist.name}{ optional name of distance }
  \item{muDatOrMod}{ a character string specifying whether
        as integration measure \eqn{mu} in Cramer-von-Mises distance,
        the empirical cdf (corresponding to argument value
        \code{"Dat"}) or the current model distribution
        (corresponding to argument value \code{"Mod"}) or a given
        integration (probability) measure / distribution \code{mu}
        (corresponding to argument value \code{"Other"}) is to be used;
         must be one of "Dat" (default) or "Mod" or "Other".
         You can specify just the initial letter; the default is \code{"Mod"}.}
  \item{mu}{ optional integration (probability) measure for CvM MDE.
             defaults to \code{NULL} and is ignored in options
             \code{muDatOrMod} in \code{"Dat"} and \code{"Mod"};
             in case \code{"Other"}, it must be of class
             \code{UnivariateDistribution}. }
  \item{paramDepDist}{logical; will computation of distance be parameter
                      dependent (see also note below)? if \code{TRUE}, distance function
                      must be able to digest a parameter \code{thetaPar}; otherwise
                      this parameter will be eliminated if present in \code{...}-argument.}
  \item{startPar}{ initial information used by \code{optimize} resp. \code{optim};
    i.e; if (total) parameter is of length 1, \code{startPar} is 
    a search interval, else it is an initial parameter value; if \code{NULL}
    slot \code{startPar} of \code{ParamFamily} is used to produce it;
    in the multivariate case, \code{startPar} may also be of class \code{Estimate},
    in which case slot \code{untransformed.estimate} is used.}
  \item{Infos}{ character: optional informations about estimator }
  \item{trafo}{ an object of class \code{MatrixorFunction} -- a transformation
  for the main parameter}
  \item{penalty}{(non-negative) numeric: penalizes non valid parameter-values}
  \item{validity.check}{logical: shall return parameter value be checked for
  validity? Defaults to yes (\code{TRUE})}
  \item{asvar.fct}{optionally: a function to determine the corresponding
    asymptotic variance; if given, \code{asvar.fct} takes arguments
    \code{L2Fam}((the parametric model as object of class \code{L2ParamFamily})) 
    and \code{param} (the parameter value as object of class 
    \code{ParamFamParameter}); arguments are called by name; \code{asvar.fct}
     may also process further arguments passed through the \code{\dots} argument}              
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{\dots}{ for the estimators: further arguments to \code{criterion} or \code{optimize}
    or \code{optim}, respectively; for \code{CvMDist2},
     these can be used e.g. by \code{E()}. }
  \item{.withEvalAsVar}{logical: shall slot \code{asVar} be evaluated
                   (if \code{asvar.fct} is given) or
                   just the call be returned?}
  \item{nmsffx}{character: a potential suffix to be appended to the estimator name.}
  \item{e1}{object of class \code{"Distribution"} or class \code{"numeric"} }
  \item{e2}{object of class \code{"Distribution"} }
  \item{.with.checkEstClassForParamFamily}{logical: Should a the end of the
     function \code{.checkEstClassForParamFamily}; defaults to \code{TRUE};
     can be switched off for computational time or because this is already
     checked in a calling wrapper function.}
}
\details{
  The argument \code{distance} has to be a (generic) function with arguments 
  the empirical data as well as an object of class \code{"Distribution"} 
  and possibly \code{\dots}; e.g. \code{KolmogorovDist} (default), 
  \code{TotalVarDist} or \code{HellingerDist}. Uses \code{\link{mceCalc}}
  for method dispatch.

  The functions \code{CvMMDEstimator}, \code{KolmogorovMDEstimator},
  \code{TotalVarMDEstimator}, and \code{HellingerMDEstimator} are aliases
  where the distance is fixed. More specifically, \code{CvMMDEstimator}
  uses Cramer-von-Mises distance, see \code{\link[distrEx]{CvMDist}}
  with integration measure \eqn{mu} either
  equal to the empirical cdf or to the current best fitting model distribution;
  the alternative is selected by argument \code{muDatOrMod}).
  As it is asymptotically linear, asymptotic variances are available.
  In case of alternative \code{"Dat"}, this variance is computed by means
  of helper function \code{.CvMMDCovarianceWithMux}, case of alternative
  \code{"Mod"} we use helper function \code{.CvMMDCovariance}. In both
  case one may use these helper function to get hand on the respective
  influence function. For covariances computed by \code{.CvMMDCovariance},
  diagnostics on the involved integrations are available if argument
  \code{diagnostic} is \code{TRUE}. Then there is attribute \code{diagnostic}
  attached to the return value, which may be inspected
  and accessed through \code{\link[distrEx:distrExIntegrate]{showDiagnostic}} and
  \code{\link[distrEx:distrExIntegrate]{getDiagnostic}}.


  \code{KolmogorovMDEstimator} uses Kolmogorov distance, see
  \code{\link[distrEx]{KolmogorovDist}}, \code{TotalVarMDEstimator},
  uses total variation distance, see \code{\link[distrEx]{TotalVarDist}}
  and \code{HellingerMDEstimator}
  uses Hellinger distance, see \code{\link[distrEx]{HellingerDist}}.

  Function \code{CvMDist2} calls \code{\link[distrEx]{CvMDist}} and
  computes the Cramer-von-Mises distance between
  distributions \code{e1} and \code{e2} with integration
  measure \code{mu} equal to \code{e2}; it is used in alternative
  \code{"Mod"} in \code{CvMMDEstimator}.
}
\value{
  The estimators return an object of S4-class \code{"MCEstimate"} which inherits from class
  \code{"Estimate"}. \code{CvMDist2} returns the respective distance.
}
\note{The distance function may be called together with a parameter \code{thetaPar}
 which is the current parameter value under consideration, i.e.; the value
 under which the model distribution is considered. Hence, if desired,
 particular distance functions could make use of this information, by, say
 computing the distance differently for different parameter values.}

\section{Theoretical Background}{It should be noted that \code{CvMMDEstimator}
  results in an asymptotically linear (hence asymptotically normal) estimator
  with an influence function which is always bounded;
  \code{HellingerMDEstimator} adapts, for growing sample size,
  the MLE estimator, hence is asymptotically efficient, while for finite
  sample size is bias robust. \code{KolmogorovMDEstimator} is square-root-n
  consistent but, due to the facetted level sets of the distance fails to
  be asymptotically normal. In the terminology of Donoho/Liu,
  \code{TotalVarMDEstimator} and \code{HellingerMDEstimator} rely on
  strong distances, while \code{CvMMDEstimator} and
  \code{KolmogorovMDEstimator} use weak distances, so the latter ensure
  protection against larger classes of contamination (simply because the
  distribution balls based on the respective distances contain more elements).
 }

\references{
    Beran, R. (1977). Minimum Hellinger distance estimates for parametric models.
    \emph{Annals of Statistics}, \bold{5}(3), 445-463.

    Donoho, D.L. and Liu, R.C. (1988). The "automatic" robustness of minimum
    distance functionals. \emph{Annals of Statistics}, \bold{16}(2), 552-586.

    Huber, P.J. (1981) \emph{Robust Statistics}. New York: Wiley.

    Parr, W.C. and Schucany, W.R. (1980). Minimum distance and robust estimation.
    \emph{Journal of the American Statistical Association}, \bold{75}(371),
    616-624.

    Rao, P.V., Schuster, E.F., and Littell, R.C. (1975).
    Estimation of Shift and Center of Symmetry Based on Kolmogorov-Smirnov
    Statistics. \emph{Annals of Statistics}, \bold{3}, 862-873.

    Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{ParamFamily-class}}, \code{\link{ParamFamily}}, 
         \code{\link{MCEstimator}}, \code{\link{MCEstimate-class}},
          \code{\link[MASS]{fitdistr}} }
\examples{
## (empirical) Data
set.seed(123)
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

## Kolmogorov(-Smirnov) minimum distance estimator
MDEstimator(x = x, ParamFamily = G, distance = KolmogorovDist)
## or
KolmogorovMDEstimator(x = x, ParamFamily = G)

## von Mises minimum distance estimator with default mu = Mod
MDEstimator(x = x, ParamFamily = G, distance = CvMDist)

\donttest{
## von Mises minimum distance estimator with default mu = Mod
MDEstimator(x = x, ParamFamily = G, distance = CvMDist,
            asvar.fct = .CvMMDCovarianceWithMux)
## or
CvMMDEstimator(x = x, ParamFamily = G)
## or
CvMMDEstimator(x = x, ParamFamily = G, muDatOrMod="Mod")

## or with data based integration measure:
CvMMDEstimator(x = x, ParamFamily = G, muDatOrMod="Dat")

## von Mises minimum distance estimator with mu = N(0,1)
MDEstimator(x = x, ParamFamily = G, distance = CvMDist, mu = Norm())
## or, with asy Var
MDEstimator(x = x, ParamFamily = G, distance = CvMDist, mu = Norm(),
            asvar.fct = function(L2Fam, param, ...){
            .CvMMDCovariance(L2Fam=L2Fam, param=param, mu=Norm(), N = 400)
            } )
## synomymous to
CvMMDEstimator(x = x, ParamFamily = G, muDatOrMod="Other", mu = Norm())

## Total variation minimum distance estimator
## gamma distributions are discretized
MDEstimator(x = x, ParamFamily = G, distance = TotalVarDist)
## or
TotalVarMDEstimator(x = x, ParamFamily = G)
## or smoothing of emprical distribution (takes some time!)
#MDEstimator(x = x, ParamFamily = G, distance = TotalVarDist, asis.smooth.discretize = "smooth")

## Hellinger minimum distance estimator
## gamma distributions are discretized
distroptions(DistrResolution = 1e-10)
MDEstimator(x = x, ParamFamily = G, distance = HellingerDist, startPar = c(1,2))
## or
HellingerMDEstimator(x = x, ParamFamily = G, startPar = c(1,2))
distroptions(DistrResolution = 1e-6) # default
## or smoothing of emprical distribution (takes some time!)
MDEstimator(x = x, ParamFamily = G, distance = HellingerDist, asis.smooth.discretize = "smooth")
}
}
\keyword{univar}
\keyword{robust}
