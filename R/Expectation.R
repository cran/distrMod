## why did I choose "L2ParamFamily and not "ParamFamily?
setMethod("E", signature(object = "L2ParamFamily", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
        return(E(object = object@distribution, fun = fun, useApply = useApply, ..., diagnostic = diagnostic))
    })
setMethod("E", signature(object = "L2ParamFamily", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
        res <- E(object = object, fun = as(fun, "EuclRandVariable"),
                 useApply = useApply, ...)
        if(diagnostic){
           diagn <- attr(res,"diagnostic")
           diagn[["call"]] <- match.call()
        }
        res <- matrix(res, nrow = nrow(fun))
        if(diagnostic){
           attr(res, "diagnostic") <- diagn
           class(attr(res,"diagnostic"))<- "DiagnosticClass"
        }
        return(res)
    })
setMethod("E", signature(object = "L2ParamFamily", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        diagn <- if(diagnostic) vector("list",nrvalues) else NULL
        for(i in 1:nrvalues){
            res[[i]] <- buf <- E(object = object, fun = fun[[i]],
                                 useApply = useApply, ..., diagnostic = diagnostic)
            if(diagnostic) diagn[[i]] <- attr(buf,"diagnostic")
        }
        if(diagnostic){
           diagn <- attr(res,"diagnostic")
           diagn[["call"]] <- match.call()
           attr(res, "diagnostic") <- diagn
           class(attr(res,"diagnostic"))<- "DiagnosticClass"
       }
        return(res)
    })
