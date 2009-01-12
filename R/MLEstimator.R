###############################################################################
## Implementation of Maximum Likelihood estimation in i.i.d. setup
###############################################################################


## Maximum-Likelihood estimator
MLEstimator <- function(x, ParamFamily, startPar = NULL, 
                        Infos, trafo = NULL, penalty = 0, ...){

    ## preparation: getting the matched call
    es.call <- match.call()
    dots <- match.call(expand.dots = FALSE)$"..."


    ## some checking
    if(!is.numeric(x))
      stop(gettext("'x' has to be a numeric vector"))   
    if(is.null(startPar)) startPar <- startPar(ParamFamily)(x,...)


    ## manipulation of the arg list to method mceCalc
    argList <- c(list(x = x, PFam = ParamFamily, startPar = startPar, 
                      penalty = penalty))
    if(missing(Infos))      Infos <- NULL
        argList <- c(argList, Infos = Infos)
    if(!is.null(dots))      argList <- c(argList, dots)

    ## call to mleCalc
    res0 <- do.call(mleCalc, argList)
    
    ## setting asymptotic variance
    asv.fct <- function(PFam = ParamFamily, param, ...) 
                        solve(FisherInfo(PFam, param = param))

    ## digesting the results of mceCalc
    res <- .process.meCalcRes(res0, PFam = ParamFamily, 
                              trafo = trafo, 
                              res.name = "Maximum likelihood estimate", 
                              call = es.call,                               
                              asvar.fct = asv.fct, 
                              ...)

    
    names(res@criterion) <- "negative log-likelihood"
    res@estimate.call <- es.call
    res@name <- "Maximum likelihood estimate"

    
    return(res)
}

 