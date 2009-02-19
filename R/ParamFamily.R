### from Matthias' thesis / ROptEst
### extended: new slot/argument modifyParam
ParamFamily <- function(name, distribution = Norm(), distrSymm,
                        modifyParam, 
                        main = main(param), nuisance = nuisance(param),
                        fixed = fixed(param), trafo = trafo(param),
                        param = ParamFamParameter(name = paste("Parameter of", 
                                      name),  main = main, nuisance = nuisance, 
                                              fixed = fixed, trafo = trafo),
                        props = character(0), startPar = NULL, makeOKPar = NULL){
    f.call <- match.call()
    if(missing(name))
        name <- "parametric family of probability measures"
    if(missing(distrSymm)) distrSymm <- NoSymmetry()
    if(missing(param)&&missing(main))
        param <- ParamFamParameter(name = "location", main = 0)
    if(missing(param)){
        argList <- list(name = paste("Parameter of", name),
                                   main = main)
        if(!missing(nuisance)) argList <- c(argList, nuisance = nuisance)                            
        if(!missing(fixed))    argList <- c(argList, fixed = fixed)                            
        if(!missing(trafo))    argList <- c(argList, trafo = trafo)                            
        param <- do.call(ParamFamParameter, argList)
        
    }
    if(missing(modifyParam))
       stop(cat(paste("Please enter a function(theta) with value a new instance of",
                      "slot distribution with corresponding parameter value theta.",
                      "example (normal location) function(theta){ Norm(mean=theta)}\n", collapse="",sep="")))
    PF <- new("ParamFamily")
    PF@name <- name
    PF@distribution <- distribution
    PF@distrSymm <- distrSymm
    PF@param <- param
    PF@props <- props
    PF@modifyParam <- modifyParam
    PF@fam.call <- f.call
    if(!is.null(startPar)) PF@startPar <- startPar
    if(!is.null(makeOKPar)) PF@makeOKPar <- makeOKPar

    return(PF)
}


## access methods
setMethod("param", "ParamFamily", function(object) object@param)
setMethod("modifyParam", "ParamFamily", 
    function(object){
        fun <- function(theta){}
        body(fun) <- substitute({ validParameter(object, param = theta); fun(theta) },
                                list(fun = object@modifyParam))
        return(fun)
    })
setMethod("startPar", "ParamFamily", 
    function(object){
        fun <- function(x,...){}
        body(fun) <- substitute(fun(x ,...),
                                list(fun = object@startPar))
        return(fun)
    })
setMethod("makeOKPar", "ParamFamily", 
    function(object){
        fun <- function(x,...){}
        body(fun) <- substitute(fun(x ,...),
                                list(fun = object@makeOKPar))
        return(fun)
    })
setMethod("fam.call", "ParamFamily", function(object) object@fam.call)

## wrapped access methods
setMethod("main", "ParamFamily", function(object) main(param(object)))
setMethod("nuisance", "ParamFamily", function(object) nuisance(param(object)))
setMethod("fixed", "ParamFamily", function(object) fixed(param(object)))
setMethod("trafo", signature(object = "ParamFamily", param = "missing"), 
                   function(object, param){ param0 <- object@param 
                                            return(trafo(param0))})
setMethod("trafo", signature(object = "ParamFamily", param = "ParamFamParameter"), 
   function(object, param){
        if(is.function(trafo(object))) 
             return(list(fct = trafo(object), 
                         mat = (trafo(object)(main(param)))$mat))
        else return(list(fct = function(x) trafo(object)%*%x, 
                         mat = trafo(object)))
   })  

## replace methods
#setReplaceMethod("param", "ParamFamily", 
#    function(object, value){ object@param <- value; object })
setReplaceMethod("trafo", "ParamFamily", 
    function(object, value){ 
        param <- object@param
        param@trafo <-  value
        object <- modifyModel(object, param)
        object
    })

