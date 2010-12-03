\name{distrMod-package}
\alias{distrMod-package}
\alias{distrMod}
\docType{package}

\title{
distrMod -- probability models
}
\description{
Based on the packages \pkg{distr} and \pkg{distrEx} package \pkg{distrMod} 
provides a flexible framework which allows computation of estimators like 
maximum likelihood or minimum distance estimators for probability models.
}

\details{
\tabular{ll}{
Package: \tab distrMod \cr
Version: \tab 2.3 \cr
Date: \tab 2010-12-03 \cr
Depends: \tab R(>= 2.6.0), methods, startupmsg, distr(>= 2.2), distrEx(>=
2.2), RandVar(>= 0.6.3), MASS, stats4 \cr
LazyLoad: \tab yes \cr
License: \tab LGPL-3 \cr
URL: \tab http://distr.r-forge.r-project.org/\cr
SVNRevision: \tab 699 \cr
}}

\section{Classes}{
\preformatted{
[*]: there is a generating function with the same name

##########################
ProbFamily classes
##########################
slots: [<name>(<class>)] 
name(character), distribution(Distribution),
distrSym(DistributionSymmetry), props(character)

"ProbFamily"
|>"ParamFamily"     [*] 
additional slots:
param(ParamFamParameter), modifyParam(function),
startPar(function), makeOKPar(function), fam.call(call)
|>|>"L2ParamFamily" [*]
additional slots:
L2deriv(EuclRandVarList), L2deriv.fct(function), 
L2derivSymm(FunSymmList), L2derivDistr(DistrList),
L2derivDistrSymm(DistrSymmList), FisherInfo(PosSemDefSymmMatrix),
FisherInfo.fct(function)
|>|>|>"BinomFamily" [*]
|>|>|>"PoisFamily"  [*]
|>|>|>"GammaFamily" [*]
|>|>|>"BetaFamily"  [*]
|>|>|>"L2GroupParamFamily"
additional slots:
LogDeriv(function)
|>|>|>|>"L2LocationScaleUnion"  /VIRTUAL/
additional slots:
locscalename(character)
|>|>|>|>|>"L2LocationFamily"            [*]
|>|>|>|>|>|>"NormLocationFamily"        [*] 
|>|>|>|>|>|>"GumbelLocationFamily"      [*] 
|>|>|>|>|>"L2ScaleFamily"               [*]
|>|>|>|>|>|>"NormScaleFamily"           [*] 
|>|>|>|>|>|>"ExpScaleFamily"            [*] 
|>|>|>|>|>|>"LnormScaleFamily"          [*] 
|>|>|>|>|>"L2LocationScaleFamily"       [*]
|>|>|>|>|>|>"NormLocationScaleFamily"   [*] 
|>|>|>|>|>|>"CauchyLocationScaleFamily" [*] 

##########################
ParamFamParameter
##########################
"ParamFamParameter" [*] is subclass of class "Parameter" of package "distr".
Additional slots:
main(numeric), nuisance(OptionalNumeric), fixed(OptionalNumeric),
trafo(MatrixorFunction)

##########################
Class unions
##########################
"MatrixorFunction" = union("matrix", "OptionalFunction")
"PrintDetails" = union("Estimate", "Confint", 
                   "PosSemDefSymmMatrix",
                   "ParamFamParameter", "ParamFamily")


##########################
Symmetry classes            (other classes moved to package "distr")
##########################
slots:
type(character), SymmCenter(ANY)

"Symmetry"   (from package "distr")
|>"FunctionSymmetry"
|>|>"NonSymmetric"      [*]
|>|>"EvenSymmetric"     [*]
|>|>"OddSymmetric"      [*]

list thereof
"FunSymmList"           [*]


##########################
Matrix classes              (moved to package "distr")
##########################
slots:
none
"PosSemDefSymmMatrix" [*] is subclass of class "matrix" of package "base".
|>"PosDefSymmMatrix"  [*]


##########################
Norm Classes
##########################
slots:
name(character), fct(function)

"NormType"        [*]
|>"QFNorm"        [*]
Additional slots:
QuadForm(PosSemDefSymmMatrix)
|>|>"InfoNorm"    [*]
|>|>"SelfNorm"    [*]


##########################
Bias Classes
##########################
slots:
name(character)

"BiasType"
|>"symmetricBias"   [*]
|>"onesidedBias"
Additional slots:
sign(numeric)
|>"asymmetricBias"  [*]
Additional slots:
nu(numeric)


##########################
Risk Classes
##########################
slots:
type(character)

"RiskType"
|>"asRisk"
|>|>"asCov"       [*]
|>|>"trAsCov"     [*]
|>"fiRisk"
|>|>"fiCov"       [*]
|>|>"trfiCov"     [*]
|>|>"fiHampel"    [*]
Additional slots:
bound(numeric)
|>|>"fiMSE"       [*]
|>|>"fiBias"      [*]
|>|>"fiUnOvShoot" [*]
Additional slots:
width(numeric)

Risk with Bias:
"asRiskwithBias"
slots: biastype(BiasType), normtype(NormType),
|>"asHampel"      [*]
Additional slots:
bound(numeric)
|>"asBias"        [*]
|>"asGRisk"
|>|>"asMSE"       [*]
|>|>"asUnOvShoot" [*]
Additional slots:
width(numeric)
|>|>"asSemivar"   [*]


##########################
Estimate Classes
##########################
slots:
name(character), estimate(ANY),
samplesize(numeric), asvar(OptionalMatrix),
Infos(matrix), nuis.idx(OptionalNumeric)
fixed.estimate(OptionalNumeric), 
estimate.call(call), trafo(list[of function, matrix]),
untransformed.estimate(ANY), 
untransformed.asvar(OptionalMatrix)
criterion.fct(function), method(character),

"Estimate"
|>"MCEstimate", 
Additional slots:
criterion(numeric)

##########################
Confidence interval class
##########################
slots:
type(character), confint(array),
estimate.call(call), name.estimate(character),
trafo.estimate(list[of function, matrix]), 
nuisance.estimate(OptionalNumeric)
"Confint"

}
}
\section{Methods}{
besides accessor and replacement functions, we have
methods 

\code{solve}, \code{sqrt} for matrices
\code{checkL2deriv}, \code{existsPIC} for class \code{L2ParamFamily}
\code{LogDeriv} for  class \code{L2GroupParamFamily}

\code{validParameter} for classes \code{ParamFamily}, \code{L2ScaleFamily},
\code{L2LocationFamily}, and \code{L2LocationScaleFamily}

\code{modifyModel} for the pairs of classes 
\code{L2ParamFamily} and \code{ParamFamParameter},
\code{L2LocationFamily} and \code{ParamFamParameter},
\code{L2ScaleFamily} and \code{ParamFamParameter},
\code{L2LocationScaleFamily} and \code{ParamFamParameter},
\code{GammaFamily} and \code{ParamFamParameter}, and
\code{ExpScaleFamily} and \code{ParamFamParameter}

\code{mceCalc} for the pair of classes \code{numeric} and \code{ParamFamily}

\code{mleCalc} for the pairs of classes
\code{numeric} and \code{ParamFamily},
\code{numeric} and \code{BinomFamily},
\code{numeric} and \code{PoisFamily},
\code{numeric} and \code{NormLocationFamily},
\code{numeric} and \code{NormScaleFamily}, and
\code{numeric} and \code{NormLocationScaleFamily}

\code{coerce} from class \code{MCEstimate} to class \code{mle}

\code{confint} for  class \code{Estimate}
\code{profile} for  class \code{MCEstimate}

}
\section{Functions}{

\preformatted{
Management of global options:
"distrModOptions", "distrModoptions", "getdistrModOption",

check for ker of matrix: "isKerAinKerB"

particular norms: "EuclideanNorm", "QuadFormNorm"
onesided bias: "positiveBias", "negativeBias", 

Estimators:
"Estimator", "MCEstimator", "MLEstimator", "MDEstimator"

special location/scale models:
"L2LocationUnknownScaleFamily", "L2ScaleUnknownLocationFamily"

some special normal models:
"NormScaleUnknownLocationFamily", "NormLocationUnknownScaleFamily",
}}

\note{
Some functions of packages \pkg{stats}, \pkg{base} have intentionally been masked, 
but completely retain their functionality --- see \code{distrModMASK()}.

If any of the packages \pkg{stats4}, \pkg{fBasics} is to be used  
together with \pkg{distrMod}, the latter must be attached \emph{after} any of the 
first mentioned. Otherwise \code{confint()} defined as \emph{method} 
in \pkg{distrMod} may get masked.\cr  To re-mask, you 
may use  \code{confint <- distrMod::confint}. 
See also \code{distrModMASK()}
}

\author{
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de},\cr 
Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr

\emph{Maintainer:} Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
}

\references{
M. Kohl and P. Ruckdeschel (2010):
R Package distrMod: S4 Classes and Methods for Probability Models.
Journal of Statistical Software, 35(10), 1-27. 
\url{http://www.jstatsoft.org/v35/i10/}
(see also \code{vignette("distrMod")})


P. Ruckdeschel, M. Kohl, T. Stabla, F. Camphausen (2006):
S4 Classes for Distributions, \emph{R News}, \emph{6}(2), 2-6. 
\url{http://CRAN.R-project.org/doc/Rnews/Rnews_2006-2.pdf}

A vignette for packages \pkg{distr}, \pkg{distrSim}, \pkg{distrTEst}, and 
\pkg{distrEx} is included into the mere documentation package \pkg{distrDoc} 
and may be called by 
\code{require("distrDoc");vignette("distr")}
}

\section{Start-up-Banner}{
You may suppress the start-up banner/message completely by setting 
\code{options("StartupBanner"="off")} somewhere before loading this 
package by \code{library} or \code{require} in your R-code / R-session.

If option \code{"StartupBanner"} is not defined (default) or setting    
\code{options("StartupBanner"=NULL)} or  
\code{options("StartupBanner"="complete")} the complete start-up banner is 
displayed.

For any other value of option \code{"StartupBanner"} 
(i.e., not in \code{c(NULL,"off","complete")})
only the version information is displayed.

The same can be achieved by wrapping the \code{library} or \code{require} 
call into either \code{suppressStartupMessages()} or 
\code{onlytypeStartupMessages(.,atypes="version")}. 

As for general \code{packageStartupMessage}'s, you may also suppress all
the start-up banner by wrapping the \code{library} or \code{require} 
call into \code{suppressPackageStartupMessages()} from 
\pkg{startupmsg}-version 0.5 on.
}

\section{Demos}{
Demos are available --- see \code{demo(package="distrMod")}.}
\section{Scripts}{
Example scripts are available --- see folder \file{scripts}
in the package folder to package \pkg{distrMod} in your library.
}

\section{Package versions}{
Note: The first two numbers of package versions do not necessarily reflect
 package-individual development, but rather are chosen for the
 distrXXX family as a whole in order to ease updating "depends"
 information.
}
\keyword{package}
