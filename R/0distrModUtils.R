.getLogDeriv <- function(distr){
  xs <- seq(2*q(distr)(getdistrOption("TruncQuantile")^2),
            2*q(distr)(getdistrOption("TruncQuantile")^2, lower = FALSE),
            length = getdistrOption("DefaultNrGridPoints"))
  m <- getdistrOption("DefaultNrGridPoints")%/%100
  dxs<- -d(distr)(xs, log = TRUE)
  x1 <- xs[1]; xn <- (rev(xs)[1])
  f2xs <- approxfun(x = xs, y = D2ss(xs,dxs)$y, rule = 2)
  f2x1 <- f2xs(x1); f2xn <- f2xs(xn);
  f1xs <- approxfun(x = xs, y = D1ss(xs,dxs))
  f1x1 <- f1xs(x1); f1xn <- f1xs(xn);
  f3xs <- approxfun(x = xs, y = D2ss(xs,f1xs(xs))$y, rule = 1)
  f3x1 <- median(f3xs(xs[1:m])); f3xn <- median(f3xs(rev(xs)[1:m]));
  fxs <- function(x){
       f1x0 <- f1xs(x)
       dx1 <- (x[x<x1]-x1)
       dxn <- (x[x>xn]-xn)
       f1x0[x>xn] <- f1xn+f2xn*dxn+f3xn/2*dxn^2
       f1x0[x<x1] <- f1x1+f2x1*dx1+f3x1/2*dx1^2
       return(f1x0)}
  return(fxs)
}

.show.with.sd <- function(est, s){
  ### code borrowed from print.fitdistr in  package MASS by B.D. Ripley
        digits <- getOption("digits")
        ans <- format(base::rbind(est, s), digits=digits)
        ans[1, ] <- sapply(ans[1, ], function(x) paste("", x))
        ans[2, ] <- sapply(ans[2, ], function(x) paste("(", x, ")", sep=""))
     ## only used for digits
        dn <- dimnames(ans)
        dn[[1]] <- rep("", 2)
        dn[[2]] <- paste(substring("      ", 1, 
                       (nchar(ans[2,]) - nchar(dn[[2]])) %/% 2), dn[[2]])
        dn[[2]] <- paste(dn[[2]], substring("      ", 1, 
                       (nchar(ans[2,]) - nchar(dn[[2]])) %/% 2))
        dimnames(ans) <- dn
        print(ans, quote = FALSE)
        return(invisible())
        }
 ### end of borrowed code  


.isUnitMatrix <- function(m){                                                                  
### checks whether m is unit matrix
              m.row <- nrow(m)
              isTRUE(all.equal(m, diag(m.row), check.attributes = FALSE))
              }

.validTrafo <- function(trafo, dimension){
##checks whether trafo is valid
  if(!is.function(trafo)){
    if(ncol(trafo) != dimension)
        stop("invalid transformation:\n", 
             "number of columns of 'trafo' not equal to ", 
             "dimension of the parameter")
    if(nrow(trafo) > dimension)
        stop("invalid transformation:\n",
             "number of rows of 'trafo' larger than ", 
             "dimension of the parameter")
    if(any(!is.finite(trafo)))
        stop("infinite or missing values in 'trafo'")
    }
  return(invisible())
}

##caching:
.csimpsum <- distr:::.csimpsum
### still to be tested and improved:
## covariance for minimum CvM distance estimator acc. Ri:94, pp.132-133

.CvMMDCovariance<- function(L2Fam, param, mu = distribution(L2Fam), expon=3, 
                            withplot = FALSE, withpreIC = FALSE,
                            N = getdistrOption("DefaultNrGridPoints")+1, ...){

   # preparations:
   eps <- getdistrOption("TruncQuantile")


   N1 <- 2*N+1
   odd <- (1:N1)%%2==1

   param0 <- L2Fam@param
   dim0 <- dimension(param0)
   paramP <- ParamFamParameter(name = name(param0), main = main(param),
                               trafo = diag(dim0))
   L2Fam <- modifyModel(L2Fam, paramP)

   distr <- L2Fam@distribution

   if(missing(mu)) mu <- distr

   if(is(distr,"DiscreteDistribution"))
       x.seq <-support(distr)
   else
       {if(is(distr,"AbscontDistribution")){
           x.seq0 <- seq(q(distr)(eps^expon),
                         q(distr)(eps^expon, lower = FALSE), length = N1)
           h0 <- x.seq0[1:2]%*%c(-1,1)
           x.seq <- x.seq0[odd]
          }else{ 
           x.seq <- seq(q(distr)(eps^expon),
                        q(distr)(eps^expon, lower = FALSE), length = N)
          }
       }
   if(is(mu,"DiscreteDistribution"))
       x.mu.seq <-support(distr)
   else
       {if(is(mu,"AbscontDistribution")){
           x.mu.seq0 <- seq(q(mu)(eps^expon),
                            q(mu)(eps^expon, lower = FALSE), length = N1)
           h0.mu <- x.mu.seq0[1:2]%*%c(-1,1)
           x.mu.seq <- x.mu.seq0[odd]
          }else{ 
           x.mu.seq <- seq(q(mu)(eps^expon),
                        q(mu)(eps^expon, lower = FALSE), length = N)
          }
       }

   L2deriv <- L2deriv(L2Fam)[[1]]
   ## are we working with a one-dim L2deriv or not?

   onedim <- (length(L2deriv@Map)==1)


   if(onedim){
   ## one-dim case

   ## Delta, formula (56), p. 133 [Ri:94]
   ##        Ptheta- primitive function for Lambda

   if(is(distr,"AbscontDistribution")){
      Delta0x <- sapply(x.seq0, function(x) 
                                evalRandVar(L2deriv, x)) * 
                 d(distr)(x.seq0)
      Delta0 <-  h0*.csimpsum(Delta0x)   
   }else{
      L2x  <- function(x,y)  (x<=y)*evalRandVar(L2deriv, x)
      Delta0 <- sapply(x.seq, function(Y){ fct <- function(x) L2x(x,y=Y)
                                        return(E(object=distr, fun = fct))})
   }
   Delta1 <- approxfun(x.seq, Delta0, yleft = 0, yright = 0)
   if(is(distr,"DiscreteDistribution"))         
      Delta <- function(x) Delta1(x) * (x %in% support(distr))
   else  Delta <- function(x) Delta1(x)


   ## J = Var_Ptheta Delta
   J1 <- E(object=distr, fun = Delta)
   Delta.0 <- function(x) Delta(x) - J1
   J <- E(object=distr, fun = function(x) Delta.0(x)^2)

   ### CvM-IC phi
   phi <- function(x) Delta.0(x)/J

   ## integrand phi x Ptheta in formula (51) [ibid]
   phi1 <- function(x) phi(x) * p(distr)(x)
   psi1 <- E(object = mu, fun = phi1)


   ## obtaining IC psi  (formula (51))

   if(is(mu,"AbscontDistribution")){
      phix <- function(x) phi(x)*d(mu)(x)
      psi0x <- sapply(rev(x.mu.seq0), phix)
      psi0 <-  h0.mu*rev(.csimpsum(psi0x))   
   }else{
      phixy  <- function(x,y)  (x<=y)*phi(y)
      psi0 <- sapply(x.mu.seq, function(X){ fct <- function(y) phixy(x=X,y=y)
                                        return(E(object=mu, fun = fct))})
   }
   psi.1 <- approxfun(x.mu.seq, psi0, yleft = 0, yright = rev(psi0)[1])
   if(is(distr,"DiscreteDistribution"))
         psi <- function(x) (psi.1(x)-psi1) * (x %in% support(mu))
   else  psi <- function(x) psi.1(x)-psi1

   E2 <- E(object=distr, fun = function(x) psi(x)^2)
   L2deriv <- L2Fam@L2deriv[[1]]
   ## E2 = Cov_mu (psi)

#   ### control: centering & standardization
   E1 <- E(object=distr, fun = psi )
   E3 <- E(object=distr, fun = function(x) psi(x)*evalRandVar(L2deriv, x))
   psi.0 <- function(x) psi(x) - E1
   psi.01 <- function(x) psi.0(x)/E3
   if(withplot)
       { windows()
         plot(x.seq, psi.01(x.seq),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
       }
   E4 <- E(object=distr, fun = function(x) psi.01(x)^2)
   psi.01 <- EuclRandVariable(Map = list(psi.01), Domain = Reals())

#   print(list(E2,E4,E2-E4))

      }else{

   ## multivariate case

   Dim <- length(evalRandVar(L2deriv, 1))

   ## Delta, formula (56), p. 133 [Ri:94]
   ##        Ptheta- primitive function for Lambda

   Map.Delta <- vector("list",Dim)
   
   for(i in 1:Dim)
       { if(is(distr,"AbscontDistribution")){
            fct0 <- sapply(x.seq0, function(x) L2deriv@Map[[i]](x)) * 
                           d(distr)(x.seq0)
            Delta0 <-  h0*.csimpsum(fct0)   
         }else{
            fct0 <- function(x,y) L2deriv@Map[[i]](x)*(x<=y)
            Delta0 <- sapply(x.seq, function(Y){ fct <- function(x) fct0(x,y=Y)
                                            return(E(object=distr, fun = fct))})
         }         
         Delta1 <- approxfun(x.seq, Delta0, yleft = 0, yright = 0)
         if(is(distr,"DiscreteDistribution"))
               Delta <- function(x) Delta1(x) * (x %in% support(distr))
         else  Delta <- function(x) Delta1(x)
         Map.Delta[[i]] <- Delta
         env.i <- environment(Map.Delta[[i]]) <- new.env()
         assign("i", i, env=env.i)
         assign("fct", fct, env=env.i)
         assign("fct0", fct0, env=env.i)
         assign("Delta", Delta, env=env.i)
         assign("Delta0", Delta0, env=env.i)
         assign("Delta1", Delta1, env=env.i)
         if(withplot){ 
           windows()
           plot(x.seq, sapply(x.seq,Map.Delta[[i]]),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
         }

   }
   Delta <-  EuclRandVariable(Map = Map.Delta, Domain = Reals())



   ## J = Var_Ptheta Delta
   J1 <- E(object=distr, fun = Delta)
   Delta.0 <- Delta - J1
   J <- E(object=distr, fun = Delta.0 %*%t(Delta.0))
   ### CvM-IC phi
   phi <- as(solve(J)%*%Delta.0,"EuclRandVariable")

   ## integrand phi x Ptheta in formula (51) [ibid]

   Map.phi1 <- vector("list",Dim)
   for(i in 1:Dim)
       { Map.phi1[[i]] <- function(x) evalRandVar(phi,x)[i] * p(distr)(x)
         env.i <- environment(Map.phi1[[i]]) <- new.env()
         assign("i", i, env=env.i)
         }

   phi1 <- EuclRandVariable(Map = Map.phi1, Domain = Reals())
   psi1 <- E(object=mu, fun = phi1)


   ## obtaining IC psi  (formula (51))
   Map.psi <- vector("list",Dim)
   for(i in 1:Dim)
     { if(is(mu,"AbscontDistribution")){
            fct01 <- function(x) phi@Map[[i]](x)*d(mu)(x)
            fct0 <-  sapply(rev(x.mu.seq0),fct01)
            phi0 <-  h0.mu*rev(.csimpsum(fct0))   
       }else{
            fct01 <- NULL
            fct0 <- function(x,y) evalRandVar(phi, y)[i]*(x<=y)
            phi0 <- sapply(x.mu.seq, 
                           function(X){ 
                               fct <- function(y) fct0(x = X, y)
                               return(E(object = mu, fun = fct))
                               })
       }
              
       phi0a <- approxfun(x.mu.seq, phi0, yleft = 0, yright = rev(phi0)[1])
       env.i <- environment(phi1) <- new.env()
       assign("i", i, env=env.i)
       if(is(distr,"DiscreteDistribution"))
             psi0 <- function(x) phi0a(x) * (x %in% support(mu))
       else  psi0 <- function(x) phi0a(x)

       Map.psi[[i]] <- psi0
       env.i <- environment(Map.psi[[i]]) <- new.env()
       assign("i", i, env=env.i)
       assign("fct", fct, env=env.i)
       assign("fct0", fct0, env=env.i)
       assign("psi0", psi0, env=env.i)
       assign("phi0a", phi0a, env=env.i)
       assign("phi0", phi0, env=env.i)
    }
   psi <-  EuclRandVariable(Map = Map.psi, Domain = Reals())

   E2 <- E(object=distr, fun = psi %*%t(psi))   
   ## E2 = Cov_mu (psi)

   ### control: centering & standardization
   L2deriv <- L2Fam@L2deriv[[1]]
   E1 <- E(object=distr, fun = psi )
   E3 <- E(object=distr, fun = psi %*%t(L2deriv))
   psi.0 <- psi - E1
   psi.01 <- as(solve(E3)%*%psi.0,"EuclRandVariable")
   if(withplot)
      { for(i in 1:Dim)
         { windows()
           plot(x.mu.seq, sapply(x.mu.seq,psi.01@Map[[i]]),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
         }}
   E4 <- E(object=distr, fun = psi.01 %*%t(psi.01))
   }
  E4 <- PosSemDefSymmMatrix(E4)
  
  psi <-  EuclRandVarList(psi.01)
  nms <- names(c(main(param(L2Fam)),nuisance(param(L2Fam))))
  dimnames(E4) = list(nms,nms)
  if(withpreIC) return(list(preIC=psi, Var=E4))
  else return(E4)
}

### examples:
if(FALSE){
P0 <- PoisFamily();.CvMMDCovariance(P0,par=ParamFamParameter("lambda",1), withplot=TRUE)
B0 <- BinomFamily(size=8, prob=0.3);.CvMMDCovariance(B0,par=ParamFamParameter("",.3), withplot=TRUE)
N0 <- NormLocationFamily();.CvMMDCovariance(N0,par=ParamFamParameter("",0), withplot=TRUE, N = 200)
N1 <- NormScaleFamily(); re=.CvMMDCovariance(N1,par=ParamFamParameter("",1), withICwithplot=TRUE, N = 200)
NS <- NormLocationScaleFamily();.CvMMDCovariance(NS,par=ParamFamParameter("",0:1), withplot=TRUE, N = 100)
cls <- CauchyLocationScaleFamily();.CvMMDCovariance(cls,par=ParamFamParameter("",0:1), withplot=TRUE, N = 200)
Els <- L2LocationScaleFamily(loc = 0, scale = 1,
                  name = "Laplace Location and scale family",
                  centraldistribution = DExp(),
                  LogDeriv = function(x)  sign(x),
                  FisherInfo = diag(2),
                  trafo = trafo)
.CvMMDCovariance(Els,par=ParamFamParameter("",0:1), withplot=TRUE, N = 100)
}

