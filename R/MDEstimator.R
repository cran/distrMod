###############################################################################
## Implementation of minimum distance estimation
###############################################################################
MDEstimator <- function(x, ParamFamily, distance = KolmogorovDist,
                        dist.name,  paramDepDist = FALSE,
                        startPar = NULL,  Infos, 
                        trafo = NULL, penalty = 1e20,
                        validity.check = TRUE, asvar.fct, na.rm = TRUE,
                        ..., .withEvalAsVar = TRUE){

    ## preparation: getting the matched call
    es.call <- match.call()
    dots <- match.call(expand.dots = FALSE)$"..."

    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)

    ## some checking
    if(!is.numeric(x))
      stop(gettext("'x' has to be a numeric vector"))   
    if(is.null(startPar)) startPar <- startPar(ParamFamily)(x,...)
    if(missing(dist.name))
      dist.name <- names(distance(x, ParamFamily@distribution))

    if(paramDepDist) dots$thetaPar <-NULL

    distanceFctWithoutVal <- function(e1,e2,check.validity=NULL,...) distance(e1,e2,...)
    ## manipulation of the arg list to method mceCalc
    argList <- c(list(x = x, PFam = ParamFamily, criterion = distanceFctWithoutVal,
                   startPar = startPar, penalty = penalty, 
                   crit.name = dist.name, withthetaPar = paramDepDist))

    if(missing(validity.check)) validity.check <- TRUE
       argList$validity.check <- validity.check
    if(missing(Infos))      Infos <- NULL
    argList <- c(argList, Infos = Infos, check.validity = validity.check )
    if(!is.null(dots))      argList <- c(argList, dots)
    ## call to mceCalc
    res0 <- do.call(mceCalc, argList)
    ## digesting the results of mceCalc
    names(res0$criterion) <- dist.name

    argList <- c(list(res0, PFam = ParamFamily, 
                              trafo = trafo, 
                              res.name = paste("Minimum", dist.name, 
                                               "estimate", sep = " "), 
                              call = quote(es.call),
                              .withEvalAsVar = .withEvalAsVar,
                              check.validity = validity.check))

    if(!missing(asvar.fct))   argList <- c(argList, asvar.fct = asvar.fct)
    if(!is.null(dots))  argList <- c(argList, dots)
    if(!validity.check %in% names(argList))
       argList$validity.check <- TRUE

    ## digesting the results of mceCalc
    res <- do.call(.process.meCalcRes, argList)

    res@completecases <- completecases
    return(.checkEstClassForParamFamily(ParamFamily,res))
}

CvMMDEstimator <- function(x, ParamFamily, paramDepDist = FALSE,
                           startPar = NULL, Infos,
                           trafo = NULL, penalty = 1e20,
                           validity.check = TRUE, asvar.fct = .CvMMDCovariance, 
                           na.rm = TRUE, ..., .withEvalAsVar = TRUE){
  MDEstimator(x = x, ParamFamily = ParamFamily, distance = CvMDist,
              paramDepDist = paramDepDist, startPar = startPar,  Infos = Infos,
              trafo = trafo, penalty = penalty, validity.check = validity.check,
              asvar.fct = asvar.fct, na.rm = na.rm,
              ..., .withEvalAsVar = .withEvalAsVar)
}

KolmogorovMDEstimator <- function(x, ParamFamily, paramDepDist = FALSE,
                           startPar = NULL, Infos,
                           trafo = NULL, penalty = 1e20,
                           validity.check = TRUE, asvar.fct, na.rm = TRUE, ...,
                           .withEvalAsVar = TRUE){
  MDEstimator(x = x, ParamFamily = ParamFamily, distance = KolmogorovDist,
              paramDepDist = paramDepDist, startPar = startPar,  Infos = Infos,
              trafo = trafo, penalty = penalty, validity.check = validity.check,
              asvar.fct = asvar.fct, na.rm = na.rm,
              ..., .withEvalAsVar = .withEvalAsVar)
}

TotalVarMDEstimator <- function(x, ParamFamily, paramDepDist = FALSE,
                           startPar = NULL, Infos,
                           trafo = NULL, penalty = 1e20,
                           validity.check = TRUE, asvar.fct, na.rm = TRUE, ...,
                           .withEvalAsVar = TRUE){
  MDEstimator(x = x, ParamFamily = ParamFamily, distance = TotalVarDist,
              paramDepDist = paramDepDist, startPar = startPar,  Infos = Infos,
              trafo = trafo, penalty = penalty, validity.check = validity.check,
              asvar.fct = asvar.fct, na.rm = na.rm,
              ..., .withEvalAsVar = .withEvalAsVar)
}

HellingerMDEstimator <- function(x, ParamFamily, paramDepDist = FALSE,
                           startPar = NULL, Infos,
                           trafo = NULL, penalty = 1e20,
                           validity.check = TRUE, asvar.fct, na.rm = TRUE, ...,
                           .withEvalAsVar = TRUE){
  MDEstimator(x = x, ParamFamily = ParamFamily, distance = HellingerDist,
              paramDepDist = paramDepDist, startPar = startPar,  Infos = Infos,
              trafo = trafo, penalty = penalty, validity.check = validity.check,
              asvar.fct = asvar.fct, na.rm = na.rm,
              ..., .withEvalAsVar = .withEvalAsVar)
}

