###############################################################################
## Implementation of minimum distance estimation
###############################################################################
MDEstimator <- function(x, ParamFamily, distance = KolmogorovDist,
                        dist.name,  paramDepDist = FALSE,
                        startPar = NULL,  Infos, 
                        trafo = NULL, penalty = 1e20,
                        validity.check = TRUE, asvar.fct, na.rm = TRUE,
                        ..., .withEvalAsVar = TRUE, nmsffx = "",
                        .with.checkEstClassForParamFamily = TRUE){

    ## preparation: getting the matched call
    es.call <- match.call()
    dots <- match.call(expand.dots = FALSE)$"..."

    #distfc <- paste(deparse(substitute(distance)))

    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)

    ## some checking
    if(!is.numeric(x))
      stop(gettext("'x' has to be a numeric vector"))   
    if(is.null(startPar)) startPar <- startPar(ParamFamily)(x,...)

    isCvM <- FALSE
    if(missing(dist.name)){
       dist.name0 <- names(distance(x, ParamFamily@distribution))
#       print(dist.name0)
#       print(str(dist.name0))
       dist.name <- gsub("(.+distance).+","\\1", dist.name0)
       nmsffx <- paste(
           gsub(".+distance","",gsub("(.+distance) (.+)","\\2", dist.name0)),
           nmsffx, collapse=" ")
       if(isTRUE(all.equal(distance, CvMDist2))){
          dist.name <- "CvM distance"
          nmsffx <- paste("( mu = model distr. )",nmsffx, collapse=" ")
          isCvM <- TRUE
       }
       if(isTRUE(all.equal(distance,CvMDist))&&is.null(dots$mu)){
          dist.name <- "CvM distance"
          nmsffx <- paste("( mu = emp. cdf )",nmsffx, collapse=" ")
          isCvM <- TRUE
       }
       if(isTRUE(all.equal(distance,CvMDist))&&!is.null(dots$mu)){
          muc <- paste(deparse((dots$mu)))
          dots$mu <- eval(dots$mu)
          dist.name <- "CvM distance"
          nmsffx <- paste("( mu = ", muc, ")", nmsffx, collapse=" ")
          isCvM <- TRUE
       }
    }

    toClass <- "MDEstimate"
    if(any(grepl("CvMDist", paste(deparse(substitute(distance)))))) isCvM <- TRUE

    if(isCvM) toClass <- "CvMMDEstimate"

    if(paramDepDist) dots$thetaPar <-NULL

    distanceFctWithoutVal <- function(e1,e2,check.validity=NULL,...)
                     distance(e1,e2,...)
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
    argList <- c(argList, x = x)
    if(any(nmsffx!="")) argList <- c(argList, nmsffx = nmsffx)
    argList$toClass <- toClass

    ## digesting the results of mceCalc
    res <- do.call(.process.meCalcRes, argList)

    res@completecases <- completecases
    if(.with.checkEstClassForParamFamily)
         res <- .checkEstClassForParamFamily(ParamFamily,res)
    return(res)
}

CvMMDEstimator <- function(x, ParamFamily, muDatOrMod = c("Mod", "Dat", "Other"),
                           mu = NULL,
                           paramDepDist = FALSE,
                           startPar = NULL, Infos,
                           trafo = NULL, penalty = 1e20,
                           validity.check = TRUE, asvar.fct = .CvMMDCovariance, 
                           na.rm = TRUE, ..., .withEvalAsVar = TRUE,
                           nmsffx = "", .with.checkEstClassForParamFamily = TRUE){

  muDatOrMod <- match.arg(muDatOrMod)
  if(muDatOrMod=="Dat") {
     CvMDist0 <- CvMDist
     estnsffx <- "( mu = emp. cdf )"
     if(missing(asvar.fct)) asvar.fct <- .CvMMDCovarianceWithMux
  }else{
     if(muDatOrMod=="Mod") {
        CvMDist0 <- CvMDist2
        estnsffx <- "( mu = model distr. )"
        if(missing(asvar.fct)) asvar.fct <- .CvMMDCovariance
     }else{
        if(missing(mu)||is.null(mu))
           stop(gettextf("This choice of 'muDatOrMod' requires a non-null 'mu'"))
        muc <- paste(deparse(substitute(mu)))
        CvMDist0 <- function(e1,e2,... ) CvMDist(e1, e2, mu = mu, ...)
        estnsffx <- paste("( mu = ", muc, ")")
        if(missing(asvar.fct))
            asvar.fct <- function(L2Fam, param, N = 400, rel.tol=.Machine$double.eps^0.3,
                            TruncQuantile = getdistrOption("TruncQuantile"),
                            IQR.fac = 15, ...){
               .CvMMDCovariance(L2Fam=L2Fam, param=param, mu=eval(mu),
                                withplot = FALSE, withpreIC = FALSE,
                                N = N, rel.tol=rel.tol, TruncQuantile = TruncQuantile,
                                IQR.fac = IQR.fac, ...)}
     }
  }

  res <- MDEstimator(x = x, ParamFamily = ParamFamily, distance = CvMDist0,
              paramDepDist = paramDepDist, startPar = startPar,  Infos = Infos,
              trafo = trafo, penalty = penalty, validity.check = validity.check,
              asvar.fct = asvar.fct, na.rm = na.rm,
              ..., .withEvalAsVar = .withEvalAsVar,
              .with.checkEstClassForParamFamily = FALSE)
#  print(list(estnsffx, nmsffx))
  res@name <- paste("Minimum CvM distance estimate", estnsffx, nmsffx, collapse="")
  res@estimate.call <- match.call()
  if(.with.checkEstClassForParamFamily)
         res <- .checkEstClassForParamFamily(ParamFamily,res)
  return(res)
}

KolmogorovMDEstimator <- function(x, ParamFamily, paramDepDist = FALSE,
                           startPar = NULL, Infos,
                           trafo = NULL, penalty = 1e20,
                           validity.check = TRUE, asvar.fct, na.rm = TRUE, ...,
                           .withEvalAsVar = TRUE, nmsffx = "",
                           .with.checkEstClassForParamFamily = TRUE){
  res <- MDEstimator(x = x, ParamFamily = ParamFamily, distance = KolmogorovDist,
              paramDepDist = paramDepDist, startPar = startPar,  Infos = Infos,
              trafo = trafo, penalty = penalty, validity.check = validity.check,
              asvar.fct = asvar.fct, na.rm = na.rm,
              ..., .withEvalAsVar = .withEvalAsVar, nmsffx = nmsffx,
              .with.checkEstClassForParamFamily = FALSE)
  res@estimate.call <- match.call()
  if(.with.checkEstClassForParamFamily)
         res <- .checkEstClassForParamFamily(ParamFamily,res)
  return(res)
}

TotalVarMDEstimator <- function(x, ParamFamily, paramDepDist = FALSE,
                           startPar = NULL, Infos,
                           trafo = NULL, penalty = 1e20,
                           validity.check = TRUE, asvar.fct, na.rm = TRUE, ...,
                           .withEvalAsVar = TRUE, nmsffx = "",
                           .with.checkEstClassForParamFamily = TRUE){
  res <- MDEstimator(x = x, ParamFamily = ParamFamily, distance = TotalVarDist,
              paramDepDist = paramDepDist, startPar = startPar,  Infos = Infos,
              trafo = trafo, penalty = penalty, validity.check = validity.check,
              asvar.fct = asvar.fct, na.rm = na.rm,
              ..., .withEvalAsVar = .withEvalAsVar, nmsffx = nmsffx,
              .with.checkEstClassForParamFamily = FALSE)
  res@estimate.call <- match.call()
  if(.with.checkEstClassForParamFamily)
         res <- .checkEstClassForParamFamily(ParamFamily,res)
  return(res)
}

HellingerMDEstimator <- function(x, ParamFamily, paramDepDist = FALSE,
                           startPar = NULL, Infos,
                           trafo = NULL, penalty = 1e20,
                           validity.check = TRUE, asvar.fct, na.rm = TRUE, ...,
                           .withEvalAsVar = TRUE, nmsffx = "",
                           .with.checkEstClassForParamFamily = TRUE){
  res <- MDEstimator(x = x, ParamFamily = ParamFamily, distance = HellingerDist,
              paramDepDist = paramDepDist, startPar = startPar,  Infos = Infos,
              trafo = trafo, penalty = penalty, validity.check = validity.check,
              asvar.fct = asvar.fct, na.rm = na.rm,
              ..., .withEvalAsVar = .withEvalAsVar, nmsffx = nmsffx,
              .with.checkEstClassForParamFamily = FALSE)
  res@estimate.call <- match.call()
  if(.with.checkEstClassForParamFamily)
         res <- .checkEstClassForParamFamily(ParamFamily,res)
  return(res)
}

