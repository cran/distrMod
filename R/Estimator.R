###############################################################################
## Function to compute estimates
###############################################################################
Estimator <- function(x, estimator, name, Infos, asvar = NULL, nuis.idx,
                      trafo = NULL, fixed = NULL, asvar.fct, ...){

    name.est <- paste(deparse(substitute(estimator)),sep="",collapse="")     
    es.call <- match.call()
    if(missing(name))
        name <- "Some estimator"
    
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                        dimnames=list(character(0), c("method", "message")))
    else{
        Infos <- matrix(c(rep("MCEstimator", length(Infos)), Infos), ncol = 2)
        colnames(Infos) <- c("method", "message")
    }

    samplesize <- if(is.null(dim(x))) length(x) else dim(x)[2]


    estimate <- estimator(x, ...)
    
    l.e <- length(estimate)
    idx <- idm <- 1:l.e
    mat <- diag(l.e)

    name.est <- paste(name.est,idm, sep="")     

    res <- new("Estimate")

    res@samplesize <- samplesize
    res@untransformed.estimate <- estimate
    res@estimate.call <- es.call
    res@name <- name
    res@Infos <- Infos
    res@fixed <- fixed
    
    if(missing(nuis.idx)) res@nuis.idx <- NULL
    else res@nuis.idx <- nuis.idx

    if(!is.null(res@nuis.idx))
        {idx <- res@nuis.idx
         idm <- idm[-idx]
         mat <- diag(length(idm))}
    
    if(is.null(names(estimate))) names(estimate) <- name.est
    
    param <- ParamFamParameter(name = names(estimate), 
                               main = res@estimate[idm],
                               nuisance = res@estimate[idx],
                               fixed = fixed)
    
    if(missing(trafo)||is.null(trafo)) 
       res@trafo <- list(fct = function(x) 
                               list(fval = x, mat = mat), 
                         mat = mat)
    else {if(is.matrix(trafo))
             res@trafo <- list(fct = function(x) 
                                     list(fval = trafo %*% x, mat = trafo), 
                               mat = trafo)
          else
             res@trafo <- list(fct = trafo, mat = trafo(main(param))$mat)           
         } 

    res@estimate <- estimate[idx]
    
    asvar <- NULL
    if(!missing(asvar.fct))
       asvar <- asvar.fct(L2Fam = ParamFamily, param = param, ...)

    res@untransformed.asvar <- asvar

    if(!.isUnitMatrix(res@trafo$mat)){
       res@estimate <- res@trafo$fct(estimate)
       if(!is.null(asvar))
           res@asvar <- res@trafo$mat%*%asvar[idx,idx]%*%t(res@trafo$mat)
    }


    return(res)
}
