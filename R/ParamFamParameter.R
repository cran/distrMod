### from Matthias' thesis / ROptEst
## generating function
ParamFamParameter <- function(name, main = numeric(0), nuisance, fixed, trafo){
    if(missing(name))
        name <- "parameter of a parametric family of probability measures"
    if(missing(nuisance))
        nuisance <- NULL
    if(missing(fixed))
        fixed <- NULL
    if(missing(trafo))
        trafo <- diag(length(main))

    dimension <- length(main) #+ length(nuisance)
    .validTrafo(trafo, dimension) ### check validity
    PFP <- new("ParamFamParameter")
    PFP@name <- name
    PFP@main <- main
    PFP@nuisance <- nuisance
    PFP@fixed <- fixed
    PFP@trafo <- trafo

    return(PFP)
}

## access methods
setMethod("main", "ParamFamParameter", function(object) object@main)
setMethod("nuisance", "ParamFamParameter", function(object) object@nuisance)
setMethod("fixed", "ParamFamParameter", function(object) object@fixed)
setMethod("trafo", signature(object = "ParamFamParameter", param = "missing"),
 function(object, param){ 
   if(is.function(object@trafo)) {
        main0 <- main(object)
        retv <- object@trafo(main0)
        return(retv$mat)}
   else return(object@trafo)})  

## replace methods
setReplaceMethod("main", "ParamFamParameter", 
    function(object, value){ 
        object@main <- value
        dimension <- length(object@main) # + length(object@nuisance)
        .validTrafo(object@trafo, dimension)
        object
    })
setReplaceMethod("nuisance", "ParamFamParameter", 
    function(object, value){ 
        object@nuisance <- value
        object
    })
setReplaceMethod("fixed", "ParamFamParameter", 
    function(object, value){ 
        object@fixed <- value
        object
    })
setReplaceMethod("trafo", "ParamFamParameter", 
    function(object, value){ 
        dimension <- length(object@main)# + length(object@nuisance)
        .validTrafo(value, dimension) ### check validity
        object@trafo <- value
        object
    })

## method length
setMethod("length", "ParamFamParameter", 
    function(x){ length(x@main) + length(x@nuisance) })

## method dimension
setMethod("dimension", "ParamFamParameter", function(object) length(object@main))
