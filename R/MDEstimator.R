###############################################################################
## Implementation of minimum distance estimation
###############################################################################
MDEstimator <- function(x, ParamFamily, distance = KolmogorovDist, dist.name, 
                        startPar = NULL,  Infos, 
                        trafo = NULL, penalty = 0, asvar.fct, ...){

    ## preparation: getting the matched call
    es.call <- match.call()
    dots <- match.call(expand.dots = FALSE)$"..."


    ## some checking
    if(!is.numeric(x))
      stop(gettext("'x' has to be a numeric vector"))   
    if(is.null(startPar)) startPar <- startPar(ParamFamily)(x,...)
    if(missing(dist.name))
      dist.name <- names(distance(x, ParamFamily@distribution))



    ## manipulation of the arg list to method mceCalc
    argList <- c(list(x = x, PFam = ParamFamily, criterion = distance, 
                   startPar = startPar, penalty = penalty, 
                   crit.name = dist.name))
    if(missing(Infos))      Infos <- NULL
    argList <- c(argList, Infos = Infos)
    if(!is.null(dots))      argList <- c(argList, dots)


    ## call to mceCalc
    res0 <- do.call(mceCalc, argList)

    ## digesting the results of mceCalc
    names(res0$criterion) <- dist.name
    res <- .process.meCalcRes(res0, PFam = ParamFamily, 
                              trafo = trafo, 
                              res.name = paste("Minimum", dist.name, 
                                               "estimate", sep = " "), 
                              call = es.call, asvar.fct = asvar.fct, ...)

    return(res)
}

