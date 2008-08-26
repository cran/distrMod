################################################################################
################################################################################
#################### mle&mceCalc methods ####################################
################################################################################
################################################################################

##helper functions
setMethod("samplesize", signature="numeric", function(object){
    if(is.null(dim(object)))
         return(length(object))
    else return(dim(object)[2])
    })

meRes <- function(x, estimate, criterion.value, param, crit.fct,
                  method = "explicit solution",
                  crit.name = "Maximum Likelihood", Infos)
        return(list(estimate = estimate, criterion = criterion.value, 
                    param = param, crit.fct = crit.fct, method = method, 
                    crit.name = crit.name, Infos = Infos, 
                    samplesize = samplesize(x)))

get.criterion.fct <- function(theta, Data, ParamFam, criterion, fun, ...){

    ### function to produce a function criterion.fct for profiling /
    ##  filling slot 'minuslogl' in object coerced to class mle:
    ##
    ##  we produce a function where all coordinates of theta appear as
    ##  separate named arguments, which then calls 'fun' with these
    ##  separate arguments again stacked to one (named) vector argument;
    ##  to this end note that in S functions and lists can be coerced
    ##  into each other, i.e. as.list(function(x1=3,x2,x3,...){<body>})
    ##  becomes a list of length length(arglist)+1, where the first
    ##  components are just the named arguments, while the last is the body
    ##  coerced to a list
    ##  realized as follows:

    dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."

    l <- length(theta)
    crit.lst <- vector("list", l+1)
    crit.lst[1:l] <- theta
    names(crit.lst) <- c(names(theta),"")
    ft <-function(){
                mc <- as.list(match.call())[-1]
                th0 <- theta
                th0[names(mc)] <- mc
                th0 <- c(unlist(th0))
                do.call(fun, c(list(theta = th0, Data = Data,
                                    ParamFamily = ParamFam,
                                    criterion = criterion) ,
                        dots))
                }
    crit.lst[l+1] <- as.list(ft)[1]
    return(as.function(crit.lst))
}


################################################################################
#### default methods for mleCalc  --- uses mceCalc
################################################################################


setMethod("mleCalc", signature(x = "numeric", PFam = "ParamFamily"),
           function(x, PFam, startPar = NULL, penalty = 0, Infos  = NULL, ...){

           res <- mceCalc(x = x, PFam = PFam, 
                          criterion = .negLoglikelihood, startPar = startPar, 
                          penalty = penalty, crit.name = "neg.Loglikelihood",
                          Infos = Infos, ...)
           names(res$criterion) <- "neg.Loglikelihood"
           return(res) 
})

################################################################################
#### default method for mceCalc --- the work-horse
################################################################################

setMethod("mceCalc", signature(x = "numeric", PFam = "ParamFamily"),
           function(x, PFam, criterion, startPar = NULL, penalty = 0, 
           crit.name = "", Infos = NULL, ...){

       if(is.null(startPar)) startPar <- startPar(PFam)(x,...)

        lmx <- length(main(PFam))
        lnx <- length(nuisance(PFam))
        fixed <- fixed(PFam)
   
       fun <- function(theta, Data, ParamFamily, criterion, ...){
               vP <- validParameter(ParamFamily, theta)
               if(!vP) theta <- makeOKPar(ParamFamily)(theta)
               if(lnx)
                  names(theta) <- c(names(main(ParamFamily)),
                                    names(nuisance(ParamFamily)))
               else  names(theta) <- names(main(ParamFamily))
               crit <- criterion(Data, ParamFamily@modifyParam(theta), ...)
               critP <- crit + penalty * (1-vP)
               return(critP)}

    if(length(param(PFam)) == 1){
        res <- optimize(f = fun, interval = startPar, Data = x,
                      ParamFamily = PFam, criterion = criterion, ...)
        theta <- res$minimum
        names(theta) <- names(main(PFam))
        crit <- res$objectiv
        method <- "optimize"
    }else{
        if(is(startPar,"Estimate")) startPar <- untransformed.estimate(startPar)
        res <- optim(par = startPar, fn = fun, Data = x, ParamFamily = PFam,
                     criterion = criterion, ...)
        theta <- as.numeric(res$par)
        names(theta) <- c(names(main(PFam)),names(nuisance(PFam)))
        method <- "optim"
        crit <- res$value
    }

    idx <-      if(lnx) lmx + 1:lnx else 1:(lmx+lnx)
    nuis.idx <- if(lnx) idx else NULL
    nuis <- if(lnx) theta[-idx] else NULL
    param <- ParamFamParameter(name = names(theta), 
                               main = theta[idx],
                               nuisance = nuis,
                               fixed = fixed)    

    crit.fct <- get.criterion.fct(theta, Data = x, ParamFam = PFam, 
                                   criterion, fun, ...)
    
    return(meRes(x, theta, crit, param, crit.fct, method = method,
                 crit.name = crit.name, Infos = Infos)) 
           })

################################################################################
####### particular methods
################################################################################

setMethod("mleCalc", signature(x = "numeric", PFam = "BinomFamily"),
           function(x, PFam, ...){
           size <- size(param(distribution(PFam)))
           theta <- mean(x)/size
           ll <- -sum(dbinom(x, size=size, prob=theta, log=TRUE))
           names(ll) <- "neg.Loglikelihood"
           crit.fct <- function(prob)
                          -sum(dbinom(x, size=size(param(PFam)), prob=prob, 
                               log=TRUE))
           param <- ParamFamParameter(name = "success probability", 
                               main = c("prob" = theta),
                               fixed = c("size" = size))
           if(!hasArg(Infos)) Infos <- NULL
           return(meRes(x, theta, ll, param, crit.fct, Infos = Infos)) 
})

setMethod("mleCalc", signature(x = "numeric", PFam = "PoisFamily"),
           function(x, PFam, ...){
           theta <- mean(x)
           ll <- -sum(dpois(x, lambda=theta, log=TRUE))
           names(ll) <- "neg.Loglikelihood"
           crit.fct <- function(lambda)
                          -sum(dpois(x, lambda=lambda, log=TRUE))
           param <- ParamFamParameter(name = "lambda", 
                               main = c("lambda"=theta))
           if(!hasArg(Infos)) Infos <- NULL
           return(meRes(x, theta, ll, param, crit.fct, Infos = Infos)) 
})

setMethod("mleCalc", signature(x = "numeric", PFam = "NormLocationFamily"),
           function(x, PFam, ...){
           theta <- mean(x); sd0 <- sd(distribution(PFam))
           ll <- -sum(dnorm(x, mean=theta, sd = sd0, log=TRUE))
           names(ll) <- "neg.Loglikelihood"
           crit.fct <- function(mean)
                           -sum(dnorm(x, mean=mean, sd = sd0, log=TRUE))
           param <- ParamFamParameter(name = "location parameter", 
                               main = c("mean"=theta))
           if(!hasArg(Infos)) Infos <- NULL
           return(meRes(x, theta, ll, param, crit.fct, Infos = Infos)) 
})

setMethod("mleCalc", signature(x = "numeric", PFam = "NormScaleFamily"),
           function(x, PFam, ...){
           theta <- sd(x); mn <- mean(distribution(PFam))
           ll <- -sum(dnorm(x, mean=mn, sd = theta, log=TRUE))
           names(ll) <- "neg.Loglikelihood"
           crit.fct <- function(sd)
                         -sum(dnorm(x, mean=mn, sd = sd, log=TRUE))  
           param <- ParamFamParameter(name = "scale parameter", 
                               main = c("sd"=theta))
           if(!hasArg(Infos)) Infos <- NULL
           return(meRes(x, theta, ll, param, crit.fct, Infos = Infos))
})

setMethod("mleCalc", signature(x = "numeric", PFam = "NormLocationScaleFamily"),
           function(x, PFam, ...){
           sd0 <- sd(x); mn <- mean(x); 
           theta <- c(mn, sd0); 
           names(theta) <- c("mean", "sd")
           ll <- -sum(dnorm(x, mean = mn, sd = sd0, log = TRUE))
           names(ll) <- "neg.Loglikelihood"
           crit.fct <- function(mean,sd)
                           -sum(dnorm(x, mean=mean, sd = sd, log=TRUE))
           param <- ParamFamParameter(name = "location and scale parameter", 
                               main = theta)
           if(!hasArg(Infos)) Infos <- NULL
           return(meRes(x, theta, ll, param, crit.fct, Infos = Infos)) 
})

