\name{MCEstimator}
\alias{MCEstimator}

\title{ Function to compute minimum criterion estimates }
\description{
  The function \code{MCEstimator} provides a general way to compute 
  estimates for a given parametric family of probability measures which 
  can be obtain by minimizing a certain criterion. For instance, 
  the negative log-Likelihood in case of the maximum likelihood
  estimator or some distance between distributions like in 
  case of minimum distance estimators.
}
\usage{
MCEstimator(x, ParamFamily, criterion, crit.name, 
            startPar = NULL, Infos, trafo = NULL, 
            penalty = 1e20, validity.check = TRUE, asvar.fct, na.rm = TRUE,
            ..., .withEvalAsVar = TRUE)
}
\arguments{
  \item{x}{ (empirical) data }
  \item{ParamFamily}{ object of class \code{"ParamFamily"} }
  \item{criterion}{ function: criterion to minimize; see Details section. }
  \item{crit.name}{ optional name for criterion. }
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
  \item{\dots}{ further arguments to \code{criterion} or \code{optimize}
    or \code{optim}, respectively. }
  \item{.withEvalAsVar}{logical: shall slot \code{asVar} be evaluated
                   (if \code{asvar.fct} is given) or
                   just the call be returned?}
}
\details{
  The argument \code{criterion} has to be a function with arguments the 
  empirical data as well as an object of class \code{"Distribution"} 
  and possibly \code{\dots}. Uses \code{\link{mceCalc}}
  for method dispatch.
}
\value{
  An object of S4-class \code{"MCEstimate"} which inherits from class 
  \code{"Estimate"}.
}
\note{The criterion function may be called together with a parameter \code{thetaPar}
 which is the current parameter value under consideration, i.e.; the value
 under which the model distribution is considered. Hence, if desired,
 particular criterion functions could make use of this information, by, say
 computing the criterion differently for different parameter values.}
%\references{  }
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
        Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{ParamFamily-class}}, \code{\link{ParamFamily}}, 
         \code{\link{MCEstimate-class}} }
\examples{
## (empirical) Data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

## Maximum Likelihood estimator
## Note: you can directly use function MLEstimator!
negLoglikelihood <- function(x, Distribution){
    res <- -sum(log(Distribution@d(x)))
    names(res) <- "Negative Log-Likelihood"
    return(res)
}
MCEstimator(x = x, ParamFamily = G, criterion = negLoglikelihood)

## Kolmogorov(-Smirnov) minimum distance estimator
## Note: you can also use function MDEstimator!
MCEstimator(x = x, ParamFamily = G, criterion = KolmogorovDist, 
            crit.name = "Kolmogorov distance")

## Total variation minimum distance estimator
## Note: you can also use function MDEstimator!
## discretize Gamma distribution
MCEstimator(x = x, ParamFamily = G, criterion = TotalVarDist, 
            crit.name = "Total variation distance")

## or smooth empirical distribution (takes some time!)
#MCEstimator(x = x, ParamFamily = G, criterion = TotalVarDist, 
#            asis.smooth.discretize = "smooth", crit.name = "Total variation distance")

## Hellinger minimum distance estimator
## Note: you can also use function MDEstimator!
## discretize Gamma distribution
distroptions(DistrResolution = 1e-8)
MCEstimator(x = x, ParamFamily = G, criterion = HellingerDist, 
            crit.name = "Hellinger Distance", startPar = c(1,2))
distroptions(DistrResolution = 1e-6)

## or smooth empirical distribution (takes some time!)
#MCEstimator(x = x, ParamFamily = G, criterion = HellingerDist, 
#            asis.smooth.discretize = "smooth", crit.name = "Hellinger distance")
}
\keyword{univar}
