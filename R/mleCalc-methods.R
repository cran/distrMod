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
                  crit.name = "Maximum Likelihood", Infos, warns = "",
                  startPar = NULL, optReturn = NULL)
        return(list(estimate = estimate, criterion = criterion.value,
                    param = param, crit.fct = crit.fct, method = method, 
                    crit.name = crit.name, Infos = Infos, 
                    samplesize = samplesize(x), warns = warns,
                    startPar = startPar, optReturn = optReturn))

get.criterion.fct <- function(theta, Data, ParamFam, criterion.ff, fun, ...){

    ### function to produce a function criterion.fct for profiling /
    ##  filling slot 'minuslogl' in object coerced to class mle:
    ##
    ##  we produce a function where all coordinates of theta appear as
    ##  separate named arguments, which then calls 'fun' with these
    ##  separate arguments again stacked to one (named) vector argument;
    ##  to this end note that in S, functions and lists can be coerced
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
                                    criterion = criterion.ff) ,
                        dots))
                }
    crit.lst[l+1] <- as.list(ft)[1]
    return(as.function(crit.lst))
}


################################################################################
#### default methods for mleCalc  --- uses mceCalc
################################################################################


setMethod("mleCalc", signature(x = "numeric", PFam = "ParamFamily"),
           function(x, PFam, startPar = NULL, penalty = 1e20,
                    dropZeroDensity = TRUE, Infos  = NULL,
                    validity.check = TRUE, ...){

           if(dropZeroDensity){
              .negLoglikelihood0 <- function(x, Distribution, ...)
                          .negLoglikelihood(x, Distribution, ...,
                                            dropZeroDensity = TRUE)
           }else{
              .negLoglikelihood0 <- function(x, Distribution, ...)
                          .negLoglikelihood(x, Distribution, ...,
                                            dropZeroDensity = FALSE)
           }

           res <- mceCalc(x = x, PFam = PFam, 
                          criterion = .negLoglikelihood0, startPar = startPar,
                          penalty = penalty, crit.name = "neg.Loglikelihood",
                          Infos = Infos, validity.check = validity.check, ...)
           names(res$criterion) <- "neg.Loglikelihood"
           return(res) 
})

################################################################################
#### default method for mceCalc --- the work-horse
################################################################################

setMethod("mceCalc", signature(x = "numeric", PFam = "ParamFamily"),
           function(x, PFam, criterion, startPar = NULL, penalty = 1e20,
           crit.name = "", Infos = NULL, validity.check = TRUE,
           withthetaPar = FALSE, ...){

       mO <- NULL
       if("makeOkPar" %in% slotNames(class(PFam))) mO <- PFam@makeOKPar
       if(is.null(mO)) mO <-  function(param)param
       if(is.null(startPar)) startPar <- mO(startPar(PFam)(x,...))

        lmx <- length(main(PFam))
        lnx <- length(nuisance(PFam))
        fixed <- fixed(PFam)

## added 2018 07 30: parse the dots argument a little to filter
## out arguments which could be of interest to optim/optimize resp. to criterion
## consequence: MCEstimator can be called with (exactly) named arguments of
## optimize/optim and with (exactly) named additional arguments (from position
## 3 on) of the criterion function; in particular if check.validity is a formal.
## All other arguments are not passed on / in particular if ... is not a formal
## and if check.validity is not a formal of the criterion, it is simply ignored.

#       mceCalcDots1 <- match.call(call = sys.call(sys.parent(1)),
#                                 expand.dots = FALSE)$"..."
       mceCalcDots <- list(...)
#       cat("------------\n");print(mceCalcDots);cat("------------\n");

       filterDots <- function(dots){
          if(length(dots)){
               dotsOptIz <- NULL
               nfmlsOptiz <- NULL

               dotsNames <- names(dots)
               if(length(param(PFam)) == 1){
                  nfmlsOptIz <- names(formals(optimize))
                  nOptProh <- c("f","interval")
               }else{
                  nfmlsOptIz <- names(formals(optim))
                  nOptProh <- c("fn","par")
               }
               nfmlsOptIz <- nfmlsOptIz[! (nfmlsOptIz %in% nOptProh)]
               if(length(nfmlsOptIz))
                  dotsOptIz <- dots[dotsNames %in% nfmlsOptIz]
               if(length(dotsOptIz)==0) dotsOptIz <- NULL

               dotsRest <- dots[!(dotsNames %in% nfmlsOptIz)]
               nfmlsCrit <- names(formals(criterion))
               dotsForCrit <- NULL

               if(length(dotsRest)&&length(nfmlsCrit)>2){
                  dotsRestNames <- names(dotsRest)
                  nmprohib <- c("Data", "theta", "ParamFamily",
                                "criterionF", nfmlsOptIz)
                  nmprohib <- nmprohib[nmprohib!="..."]
                  nfmlsCrit <- nfmlsCrit[! (nfmlsCrit %in% nmprohib)]
                  dotsRestNamesProhib <- dotsRestNames %in% nmprohib

                  dotsRest <- dotsRest[!dotsRestNamesProhib]
                  if(length(dotsRest)){
                     dotsRestNames <- names(dotsRest)
                     critL <- (dotsRestNames %in% nfmlsCrit)
                     critL <- critL | ("..." %in% nfmlsCrit)
                     dotsForCrit <- dotsRest[critL]
                  }
                  if(length(dotsForCrit)==0) dotsForCrit <- NULL
               }
               dotsForOpt <- c(dotsOptIz,dotsForCrit[!names(dotsForCrit)%in% nOptProh])
               return(list(dotsForOpt=dotsForOpt, dotsCrit=dotsForCrit, dotsOnlyOpt=dotsOptIz))
          }else return(NULL)
       }

       dotsToPass <- do.call(filterDots, list(mceCalcDots))
#       print(dotsToPass)
#       print(names(dotsToPass$dotsCrit))
       allwarns <- character(0)
       fun <- function(theta, Data, ParamFamily, criterionF, ...){
               vP <- TRUE
               dotsfun <- list(...)
               names(dotsfun) <- gsub("dotsForC\\.","",names(dotsfun))
#               cat(".....\n");print(dotsfun);cat(".....\n")
#               cat("!!!!\n")
#               print(names(dotsfun))
#               print(names(dotsToPass$dotsCrit))
#               cat("!!!!\n")
               dotsForC0 <- dotsfun[names(dotsfun)%in%names(dotsToPass$dotsCrit)]
#               print(dotsForC0)
               if(validity.check) vP <- validParameter(ParamFamily, theta)
               if(is.function(penalty)) penalty <- penalty(theta)
               if(!vP) {crit0 <- penalty; theta <- mO(theta)
               }else{
                  if(lnx)
                     names(theta) <- c(names(main(ParamFamily)),
                                       names(nuisance(ParamFamily)))
                  else  names(theta) <- names(main(ParamFamily))
                  distr.new <- try(ParamFamily@modifyParam(theta), silent = TRUE)
                  argList <- list(Data, distr.new)
                  if(!is.null(dotsForC0)) argList <- c(argList, dotsForC0)
                  if(withthetaPar) argList <- c(argList, list(thetaPar = theta))
                  if(is(distr.new,"try.error")){
                      crit0 <- penalty
                      warn0 <- paste("Parameter transformation at theta = ",
                                    paste(round(theta,3),collapse=","),
                                   " threw an error;\n",  "returning starting par;\n",
                                   sep="")
                      allwarns <<- c(allwarns,warn0)
                      warning(warn0)
                  }else{crit0 <- try(do.call(what = criterionF, args = argList),
                                     silent = TRUE)
                        if(is(crit0, "try-error")){
                            crit0 <- penalty
                            warn1 <- paste("Criterion evaluation at theta = ",
                                    paste(round(theta,3),collapse=","),
                                   " threw an error;\n",  "returning starting par;\n",
                                   sep="")
                         allwarns <<- c(allwarns,warn1)
                         warning(warn1)
                         }
                  }
               }
               critP <- crit0 + penalty * (1-vP)
               return(critP)}

    if(length(param(PFam)) == 1){
        argsOptimize <- list(f = fun, interval = startPar, Data = x,
                             ParamFamily = PFam, criterionF = criterion)
        if(!is.null(dotsToPass$dotsOnlyOpt))
            argsOptimize <- c(argsOptimize, dotsToPass$dotsOnlyOpt)
        if(!is.null(dotsToPass$dotsCrit))
            argsOptimize <- c(argsOptimize, dotsForC=dotsToPass$dotsCrit)
        optres <- do.call(optimize, argsOptimize)
        theta <- optres$minimum
        names(theta) <- names(main(PFam))
        crit <- optres$objective
        method <- "optimize"
    }else{
        if(is(startPar,"Estimate")) startPar <- untransformed.estimate(startPar)
        argsOptim <- list(par = startPar, fn = fun, Data = x,
                          ParamFamily = PFam, criterionF = criterion)
        if(!is.null(dotsToPass$dotsOnlyOpt))
            argsOptim <- c(argsOptim, dotsToPass$dotsOnlyOpt)
        if(!is.null(dotsToPass$dotsCrit))
            argsOptim <- c(argsOptim, dotsForC=dotsToPass$dotsCrit)
        optres <- do.call(optim, argsOptim)
        theta <- as.numeric(optres$par)
        names(theta) <- c(names(main(PFam)),names(nuisance(PFam)))
        method <- "optim"
        crit <- optres$value
    }

    vP <- TRUE
    if(validity.check) vP <- validParameter(PFam, theta)
    if(!vP) theta <- makeOKPar(PFam)(theta)

    idx <-      if(lnx) lmx + 1:lnx else 1:(lmx+lnx)
    nuis.idx <- if(lnx) idx else NULL
    nuis <- if(lnx) theta[-idx] else NULL

    param <- .callParamFamParameter(PFam, theta, idx, nuis, fixed)

    fun2 <- function(theta, Data, ParamFamily, criterion, ...){
               dotsTP <- filterDots(list(...))
               vP <- TRUE
               if(validity.check) vP <- validParameter(ParamFamily, theta)
               if(!vP) theta <- makeOKPar(ParamFamily)(theta)
               if(lnx)
                     names(theta) <- c(names(main(ParamFamily)),
                                       names(nuisance(ParamFamily)))
               else  names(theta) <- names(main(ParamFamily))
               distr.new <- ParamFamily@modifyParam(theta)
               crit1 <- do.call(criterion, c(list(Data, distr.new),
                                dotsToPass$dotsCrit))
               return(crit1)}

    crit.fct <- get.criterion.fct(theta, Data = x, ParamFam = PFam,
                                   criterion.ff = criterion, fun2, ...)

    return(meRes(x, theta, crit, param, crit.fct, method = method,
                 crit.name = crit.name, Infos = Infos, warns= allwarns,
                 startPar = startPar, optReturn = optres))
           })


##########################################################################
# added 2018 07 30
##########################################################################
if(FALSE){
require(distrMod)
nF <- NormLocationScaleFamily()
negLoglikelihood <- function(x, Distribution){
  res <- -sum(log(Distribution@d(x)))
  names(res) <- "Negative Log-Likelihood"
  return(res)
}
negLoglikelihood2 <- function(x, Distribution, ...){
  dots <- list(...)
  print(dots)
  res <- -sum(log(Distribution@d(x)))
  names(res) <- "Negative Log-Likelihood"
  return(res)
}
negLoglikelihood3 <- function(x, Distribution, check.validity, fn=3){
  print(c(chk=check.validity))
  print(c(fn=fn))
  res <- -sum(log(Distribution@d(x)))
  names(res) <- "Negative Log-Likelihood"
  return(res)
}
set.seed(123)
x <- rnorm(10)
MCEstimator(x = x, ParamFamily = nF, criterion = negLoglikelihood)
re <- MCEstimator(x = x, ParamFamily = nF, criterion = negLoglikelihood,
            hessian = TRUE)
re
optimReturn(re)

MCEstimator(x = x, ParamFamily = nF, criterion = negLoglikelihood2,
            fups="fu")
fo <- list(a="fu",c=list(b="e",3))
refo <- MCEstimator(x = x, ParamFamily = nF, criterion = negLoglikelihood2,
            fups=fo)
optimReturn(refo)

re2 <- MCEstimator(x = x, ParamFamily = nF, criterion = negLoglikelihood2,
            fups="fu", hessian = TRUE, fn="LU")
re2
optimReturn(re2)
MCEstimator(x = x, ParamFamily = nF, criterion = negLoglikelihood3)
MCEstimator(x = x, ParamFamily = nF, criterion = negLoglikelihood3, fn="LU")

## this shows how to do validity checks every fourth evaluation
count <- 0
negLoglikelihood4 <- function(x, Distribution, check.validity){
  count <<- count +1
  if(count %% 4==0)print(c(chk=check.validity))
  print(count)
  res <- -sum(log(Distribution@d(x)))
  names(res) <- "Negative Log-Likelihood"
  return(res)
}
MCEstimator(x = x, ParamFamily = nF, criterion = negLoglikelihood4)

}
##########################################################################
# end added 2018 07 30
##########################################################################


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
           n <- length(x)
           theta <- sqrt((n-1)/n)*sd(x); mn <- mean(distribution(PFam))
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
           n <- length(x)
           sd0 <- sqrt((n-1)/n)*sd(x); mn <- mean(x); 
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

