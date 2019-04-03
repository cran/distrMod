.CvMMDCovarianceWithMux <- function(L2Fam, param, withplot = FALSE, withpreIC = FALSE,
                            N = 1021, rel.tol=.Machine$double.eps^0.3,
                            TruncQuantile = getdistrOption("TruncQuantile"),
                            IQR.fac = 15, ..., x=NULL){
   mu <- distribution(L2Fam)
   if(!is.null(x)) mu <- DiscreteDistribution(x)
   .CvMMDCovariance(L2Fam=L2Fam, param=param, mu=mu,
                    withplot = withplot, withpreIC = withpreIC,
                    N = N, rel.tol=rel.tol, TruncQuantile = TruncQuantile,
                    IQR.fac = IQR.fac, ...)
}

CvMDist2 <- function(e1,e2,... ) {res <- CvMDist(e1, e2, mu = e2, ...)
  e2c <- paste(deparse(substitute(e2)))
  if(length(e2c) == 1)
     names(res) <- paste("CvM distance ( mu =", e2c,")")
  return(res)}

### 20180805: new function to compute asCov of CvM-MDE
#             which for the primitive functions uses integration on [0,1]
#             via quantile transformation
.CvMMDCovariance<- function(L2Fam, param, mu = distribution(L2Fam),
                            withplot = FALSE, withpreIC = FALSE,
                            N = 1021, rel.tol=.Machine$double.eps^0.3,
                            TruncQuantile = getdistrOption("TruncQuantile"),
                            IQR.fac = 15,
                            ..., diagnostic = FALSE){
   # preparations:
   dots <- list(...)

   dotsInt <- .filterEargs(dots)
   dotsInt[["f"]] <- NULL
   dotsInt[["lower"]] <- NULL
   dotsInt[["upper"]] <- NULL
   dotsInt[["stop.on.error"]] <- NULL
   dotsInt[["distr"]] <- NULL
   dotsInt[["diagnostic"]] <- NULL
   dotsInt[["useApply"]] <- NULL

   .useApply <- FALSE
   if(!is.null(dots$useApply)) .useApply <- dots$useApply

   if(missing(TruncQuantile)||TruncQuantile>1e-7) TruncQuantile <- 1e-8

   N.1 <- round(0.2*N)
   N.3 <- N.1
   N.2 <- N-N.1-N.3

   N1 <- 2*N+1
   N1.1 <- 2*N.1+1
   N1.2 <- 2*N.2+1
   N1.3 <- 2*N.3+1
   odd   <- (1:N1)%%2==1
   odd.1 <- (1:N1.1)%%2==1
   odd.2 <- (1:N1.2)%%2==1
   odd.3 <- (1:N1.3)%%2==1

   param0 <- L2Fam@param
   dim0 <- dimension(param0)

   paramP <- param0
   paramP@main <- main(param)
   paramP@trafo <- diag(dim0)
   L2Fam <- modifyModel(L2Fam, paramP)

   diagn <- if(diagnostic) list(call=match.call()) else NULL
   distr <- L2Fam@distribution

   ### get a sensible integration range:
   low <- TruncQuantile
   up <- 1-TruncQuantile


   if(is(distr,"DiscreteDistribution")){
       x.seq <-support(distr)
       if(length(x.seq) > 1000){}
       prob <- d(distr)(x.seq)
   }else{
        if(is(distr,"AbscontDistribution")){

           ## split up the integration range into
           # .1 = lower tail, .2 mid range, .3 upper tail

           x.seq0 <- seq(0, 1, length = N1)
           h0 <- diff(x.seq0[1:2])
           x.seq0.1 <- seq(0, 1, length = N1.1)
           h0.1 <- diff(x.seq0.1[1:2])
           x.seq0.2 <- seq(0, 1, length = N1.2)
           h0.2 <- diff(x.seq0.2[1:2])
           x.seq <- x.seq0[odd]
           x.seq.1 <- low+(1-low)*x.seq0.1/100
           x.seq.3 <- 1-rev(x.seq.1)
           x.seq.la <- rev(x.seq.1)[1]
           del <- 1-2*(x.seq.la+(h0.1/100+h0.2)/2)
           x.seq.2l <- x.seq.la+(h0.1/100+h0.2)/2+del*x.seq0.2
           x.seq.2r <- 1-rev(x.seq.2l)
           x.seq.2 <- (x.seq.2l+x.seq.2r)/2
           x.seq.a <- c(x.seq.1[odd.1],x.seq.2[odd.2],x.seq.3[odd.3])
          }else{
           x.seq <- seq(low,up, length = N)
          }
       }
   if(is(mu,"DiscreteDistribution")){
       x.mu.seq <- support(mu)
       if(length(x.mu.seq) > 1000){}
       prob.mu <- d(mu)(x.mu.seq)
   }else{
        if(is(mu,"AbscontDistribution")){


           x.mu.seq0 <- seq(0, 1, length = N1)
           h0.mu <- diff(x.mu.seq0[1:2])
           x.mu.seq0.1 <- seq(0, 1, length = N1.1)
           h0.mu.1 <- diff(x.mu.seq0.1[1:2])
           x.mu.seq0.2 <- seq(0, 1, length = N1.2)
           h0.mu.2 <- diff(x.mu.seq0.2[1:2])
           x.mu.seq <- x.mu.seq0[odd]
           x.mu.seq.1 <- low+(1-low)*x.mu.seq0.1/100
           x.mu.seq.3 <- 1-rev(x.mu.seq.1)
           x.mu.seq.la <- rev(x.mu.seq.1)[1]
           del.mu <- 1-2*(x.mu.seq.la+(h0.mu.1/100+h0.mu.2)/2)
           x.mu.seq.2l <- x.mu.seq.la+(h0.mu.1/100+h0.mu.2)/2+del.mu*x.mu.seq0.2
           x.mu.seq.2r <- 1-rev(x.mu.seq.2l)
           x.mu.seq.2 <- (x.mu.seq.2l+x.mu.seq.2r)/2
           x.mu.seq.a <- c(x.mu.seq.1[odd.1],x.mu.seq.2[odd.2],x.mu.seq.3[odd.3])
           x.mu.seq.b <- c(x.mu.seq.1,x.mu.seq.2,x.mu.seq.3)
            iN.mu.1 <- seq(N1.1)
            iN.mu.2 <- N1.1+seq(N1.2)
            iN.mu.3 <- N1.1+N1.2+seq(N1.3)
          }else{
           x.mu.seq <- seq(low, up, length = N)
          }
   }

   L2deriv.0 <- L2deriv(L2Fam)[[1]]
#   y.seq <- sapply(x.seq, function(x) evalRandVar(L2deriv, x))
#   plot(x.seq[!is.na(y.seq)],y.seq ,type="l")

   ## are we working with a one-dim L2deriv or not?

   onedim <- (length(L2deriv.0@Map)==1)


   myint <- function(f,...,diagnostic0 = FALSE){
      dotsFun <- .filterFunargs(dots,f)
      fwD <- function(x) do.call(f, c(list(x),dotsFun))
      do.call(distrExIntegrate, c(list(f=fwD, lower=0, upper=1,
                       stop.on.error=FALSE, distr=Unif(), diagnostic=diagnostic0,
                       dotsInt)))
   }

   if(onedim){
   ## one-dim case

   ## Delta, formula (56), p. 133 [Ri:94]
   ##        Ptheta- primitive function for Lambda

   if(is(distr,"AbscontDistribution")){
      fqx <- function(x){qx <- q.l(distr)(x)
                         evalRandVar(L2deriv.0, as.matrix(qx))[,,1]
                        }
      Delta0x.1 <- fqx(x.seq.1)
      Delta0x.2 <- fqx(x.seq.2)
      Delta0x.3 <- fqx(x.seq.3)
      Delta0.1 <-  h0/100*.csimpsum(Delta0x.1)
      Delta0.2 <-  rev(Delta0.1)[1]+h0*.csimpsum(Delta0x.2)
      Delta0.3 <-  rev(Delta0.2)[1]+h0/100*.csimpsum(Delta0x.3)
      Delta0 <- c(Delta0.1,Delta0.2,Delta0.3)
      Delta1.q <- approxfun(x.seq.a, Delta0, yleft = 0, yright = 0)
      J1 <- myint(f=Delta1.q, diagnostic0=diagnostic)
      if(diagnostic) diagn$J1 <- attr(J1,"diagnostic")
      Delta.0 <- function(x) Delta1.q(p(distr)(x))-J1
      J <- myint(f=function(x) (Delta1.q(x)-J1)^2, diagnostic0=diagnostic)
      if(diagnostic) diagn$J <- attr(J,"diagnostic")
  }else{
      if(is(distr,"DiscreteDistribution")){
         L2x <- evalRandVar(L2deriv.0, as.matrix(x.seq))[,,1]
         L2xdx <- L2x*prob
         Delta0 <- cumsum(L2xdx)
         J1 <- sum(Delta0*prob)
         Delta <- Delta0-J1
         J <- sum((Delta)^2*prob)
         Delta.0 <- approxfun(x.seq, Delta, yleft = 0, yright = 0)
         Delta <- Delta/J
      }else{
         L2x  <- function(x,y)  (x<=y)*evalRandVar(L2deriv.0, as.matrix(x))[,,1]
         diagn0 <- list()
         xseq.i <- quantile(seq(x.seq),c(0,0.25,.5,.75,1))
         Delta0 <- sapply(seq(x.seq), function(Y){
                     res <- do.call(E, c(list(object=distr,
                                         fun = function(x) L2x(x,y=x.seq[Y]),
                                         useApply = .useApply,
                                         diagnostic=diagnostic),dotsInt))
                     if(diagnostic) if(Y %in% xseq.i) diagn0[[paste(Y)]] <<- attr(res,"diagnostic")
                     })
         if(diagnostic) diagn$Delta0 <- diagn0
         Delta1 <- approxfun(x.seq, Delta0, yleft = 0, yright = 0)
         Delta <- Delta1
         J1 <- do.call(E, c(list(object=distr, fun = Delta, useApply = .useApply,
                            diagnostic=diagnostic), dotsInt))
         if(diagnostic) diagn$J1 <- attr(J1,"diagnostic")
         Delta.0 <- function(x) Delta(x) - J1
         J <- do.call(E, c(list(object=distr, fun = function(x) Delta.0(x)^2, useApply = .useApply,
                            diagnostic=diagnostic), dotsInt))
         if(diagnostic) diagn$J <- attr(J,"diagnostic")
      }
   }

   ### CvM-IC phi
   phi <- function(x) Delta.0(x)/J
#   print(head(sapply(x.seq,phi)))
   ## obtaining IC psi  (formula (51))

   if(is(mu,"AbscontDistribution")){
   ## integrand phi x Ptheta in formula (51) [ibid]
      phi1.q <- function(s){qs <- q.l(mu)(s)
                            return(phi(qs)*p(distr)(qs)) }
      psi1 <- myint(f=phi1.q, diagnostic0=diagnostic)
      if(diagnostic) diagn$psi1 <- attr(psi1,"diagnostic")

      phiqx <- function(x){qx <- q.l(mu)(x)
                          return(phi(qx))}
      psi0qx.1 <- phiqx(x.mu.seq.1)
      psi0qx.2 <- phiqx(x.mu.seq.2)
      psi0qx.3 <- phiqx(x.mu.seq.3)
      psi0q.3 <-  h0.mu/100*rev(.csimpsum(psi0qx.3))
      psi0q.2 <- psi0q.3[1]+h0.mu*rev(.csimpsum(psi0qx.2))
      psi0q.1 <- psi0q.2[1]+h0.mu/100*rev(.csimpsum(psi0qx.1))

      psi0q <-  c(psi0q.1,psi0q.2,psi0q.3)
      psi.q1 <- approxfun(x.mu.seq.a, psi0q, yleft = 0, yright = rev(psi0q)[1])
      psi.fct <- function(x) psi.q1(p(mu)(x))-psi1
   }else{
      if(is(mu,"DiscreteDistribution")&&is(distr,"DiscreteDistribution")){
         Delta.mu <- phi(x.mu.seq)
         pprob.mu <-  p(distr)(x.mu.seq)
         if(!all(support(mu)==support(distr))){
            L2x.mu   <- evalRandVar(L2deriv.0, as.matrix(x.mu.seq))[,,1]
         }else{
            L2x.mu <- L2x
         }
         psi1 <- sum(pprob.mu*Delta.mu*prob.mu)
         psi0 <- cumsum(rev(Delta.mu*prob.mu))
         psi0 <- rev(psi0)-psi1
         psi.fct <- approxfun(x.mu.seq, psi0, yleft = -psi1, yright = -psi1)
      }else{
   ## integrand phi x Ptheta in formula (51) [ibid]
         phi1 <- function(x) phi(x) * p(distr)(x)
         psi1 <- do.call(E,c(list(object = mu, fun = phi1, useApply = .useApply,
                             diagnostic = diagnostic),dotsInt))
         if(diagnostic) diagn$psi1 <- attr(psi1,"diagnostic")

         phixy  <- function(x,y)  (x<=y)*phi(y)
         diagn1 <- list()
         x.mu.seq.i <- quantile(seq(x.mu.seq),c(0,0.25,.5,.75,1))
         psi0 <- sapply(seq(x.mu.seq), function(X){
                   fct <- function(y) phixy(x=x.mu.seq[X],y=y)
                   res <- do.call(E,c(list(object=mu, fun = fct,
                                    useApply = .useApply), dotsInt))
                   if(diagnostic) if(X %in% x.mu.seq.i)
                       diagn1[[paste(X)]] <<- attr(res,"diagnostic")
                   return(res)})

         psi.1 <- approxfun(x.mu.seq, psi0, yleft = 0, yright = rev(psi0)[1])
         psi.fct <- function(x) psi.1(x)-psi1
      }
   }
 #  print(psi0)
   if(is(distr,"AbscontDistribution")){
      psi.q <- function(x){qx <- q.l(distr)(x); return(psi.fct(qx))}
   ## E2 = Cov_mu (psi)
#      E2 <- do.call(myint, c(list(f=function(x)psi.q(x)^2),dotsInt))
      E1 <- myint(f=psi.q, diagnostic0=diagnostic)
         if(diagnostic) diagn$E1 <- attr(E1,"diagnostic")
      E3 <- myint(f=function(x){
                       qx <- q.l(distr)(x)
                       L2qx <- evalRandVar(L2deriv.0,as.matrix(qx))[,,1]
                       return(psi.fct(qx)*L2qx)
                      }, diagnostic0=diagnostic)
         if(diagnostic) diagn$E3 <- attr(E3,"diagnostic")
      psi.01.f <- function(x) (psi.fct(x)-E1)/E3
      E4 <- myint(f=function(x) (psi.q(x)-E1)^2/E3^2, diagnostic0=diagnostic)
         if(diagnostic) diagn$E4 <- attr(E4,"diagnostic")
  }else{
      if(is(distr,"DiscreteDistribution")){
   ## E2 = Cov_mu (psi)
#         E2 <- sum(psi0^2*prob)
         psi0 <-  psi.fct(x.seq)
         E1 <- sum(psi0*prob)
         E3 <- sum(psi0*L2x*prob)
         psi.01d <- (psi0-E1)/E3
         E4 <- sum(psi.01d^2*prob)
         psi.01.f <- function(x) (psi.fct(x)-E1)/E3*liesInSupport(distr,x)
      }else{
   ## E2 = Cov_mu (psi)
#         E2 <- E(object=distr, fun = function(x) psi(x)^2, useApply = .useApply)
         L2x  <- function(x,y)  (x<=y)*evalRandVar(L2deriv.0, as.matrix(x))[,,1]
         E1 <- do.call(E, c(list(object=distr, fun = psi.fct,
                           useApply = .useApply, diagnostic=diagnostic), dotsInt))
         if(diagnostic) diagn$E1 <- attr(E1,"diagnostic")
         E3 <- do.call(E, c(list(object=distr, fun = function(x)
                 psi.fct(x)*evalRandVar(L2deriv.0, as.matrix(x))[,,1],
                 diagnostic=diagnostic, useApply = .useApply),dotsInt))
         if(diagnostic) diagn$E3 <- attr(E3,"diagnostic")
         psi.01.f <- function(x) (psi.fct(x) - E1)/E3
         E4 <- do.call(E, c(list(object=distr, fun = function(x) psi.01.f(x)^2,
                       diagnostic=diagnostic, useApply = .useApply),dotsInt))
         if(diagnostic) diagn$E4 <- attr(E4,"diagnostic")
      }
   }

#   ### control: centering & standardization
   if(withplot)
       { devNew() #windows()
         x0.seq <- x.seq
         if(is(distr,"AbscontDistribution")) x0.seq <- q.l(distr)(x.seq)
         plot(x0.seq, psi.01.f(x0.seq),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
       }

   psi.01 <- EuclRandVariable(Map = list(psi.01.f), Domain = Reals())


      }else{

   ## multivariate case

   Dim <- length(evalRandVar(L2deriv.0, 1))

   ## Delta, formula (56), p. 133 [Ri:94]
   ##        Ptheta- primitive function for Lambda

   Map.Delta <- vector("list",Dim)
   if(diagnostic) diagn <- list()
   if(!is(distr,"AbscontDistribution") &&
      !is(distr,"DiscreteDistribution") ) diagn$Delta0 <- vector("list",Dim)
   if(!is(mu,"AbscontDistribution") &&
      !is(mu,"DiscreteDistribution") ) diagn$phi0 <- vector("list",Dim)

   for(i in 1:Dim)
       { if(is(distr,"AbscontDistribution")){
            fct.q <- function(x){qx <- q.l(distr)(x); return(L2deriv.0@Map[[i]](qx))}
            fct0.q1 <- fct.q(x.seq.1)
            fct0.q2 <- fct.q(x.seq.2)
            fct0.q3 <- fct.q(x.seq.3)
            Delta0.q1 <-  h0/100*.csimpsum(fct0.q1)
            Delta0.q2 <-  rev(Delta0.q1)[1]+h0*.csimpsum(fct0.q2)
            Delta0.q3 <-  rev(Delta0.q2)[1]+h0/100*.csimpsum(fct0.q3)
            Delta0.q <- c(Delta0.q1,Delta0.q2,Delta0.q3)
            Delta1.q <- approxfun(x.seq.a, Delta0.q, yleft = 0, yright = 0)
            Map.Delta[[i]] <- function(x) Delta1.q(p(distr)(x))
            env.i <- environment(Map.Delta[[i]]) <- new.env()
            assign("i", i, envir=env.i)
            assign("Delta1.q", Delta1.q, envir=env.i)
         }else{
            if(is(distr,"DiscreteDistribution")){
               L2x <- evalRandVar(L2deriv.0, as.matrix(x.seq))[i,,1]
               L2xdx <- L2x*prob
               Delta.0 <- cumsum(L2xdx)
               Delta.f <- approxfun(x.seq, Delta.0, yleft = 0, yright = 0)
               Map.Delta[[i]] <-  function(x) Delta.f(x)
               env.i <- environment(Map.Delta[[i]]) <- new.env()
               assign("Delta.f", Delta.f, envir=env.i)
               assign("Delta.0", Delta.0, envir=env.i)
            }else{
               fct0 <- function(x,y) L2deriv.0@Map[[i]](x)*(x<=y)
               diagn0 <- list()
               xseq.i <- quantile(seq(x.seq),c(0,0.25,.5,.75,1))
               Delta0 <- sapply(seq(x.seq),
                           function(Y){
                              fct <- function(x) fct0(x,y=x.seq[Y])
                              res <- do.call(E,c(list(object=distr, fun = fct,
                                                 useApply=.useApply),dotsInt))
                              if(diagnostic) if(Y %in% xseq.i)
                                    diagn0[[paste(Y)]] <<- attr(res,"diagnostic")
                              return(res)
                              })
               if(diagnostic) diagn[["Delta0"]][[i]] <- diagn0
               Delta1 <- approxfun(x.seq, Delta0, yleft = 0, yright = 0)
               if(is(distr,"DiscreteDistribution"))
                     Delta <- function(x) Delta1(x) * (x %in% support(distr))
               else  Delta <- function(x) Delta1(x)
               Map.Delta[[i]] <- Delta
               env.i <- environment(Map.Delta[[i]]) <- new.env()
               assign("i", i, envir=env.i)
               assign("fct", fct, envir=env.i)
               assign("fct0", fct0, envir=env.i)
               assign("Delta", Delta, envir=env.i)
               assign("Delta0", Delta0, envir=env.i)
               assign("Delta1", Delta1, envir=env.i)
             }
         }
         #print(Delta0)
         if(withplot){
           dev.new()
           x0.seq <- x.seq
           if(is(distr,"AbscontDistribution")) x0.seq <- q.l(distr)(x.seq)
           #windows()
           plot(x0.seq, sapply(x0.seq,Map.Delta[[i]]),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
         }

   }
   Delta <-  EuclRandVariable(Map = Map.Delta, Domain = Reals())



   ## J = Var_Ptheta Delta
##-t-## print(system.time({
   J1 <- E(object=distr, fun = Delta, diagnostic=diagnostic)#, useApply = .useApply)
         if(diagnostic) diagn$J1 <- attr(J1,"diagnostic")
##-t-## }))
   Delta.0 <- Delta - J1

##-t-## print(system.time({
   J <- E(object=distr, fun = Delta.0 %*%t(Delta.0), diagnostic=diagnostic)#, useApply = .useApply)
         if(diagnostic) diagn$J <- attr(J,"diagnostic")
##-t-## }))
   ### CvM-IC phi
   phi <- as(distr::solve(J)%*%Delta.0,"EuclRandVariable")

   ## integrand phi x Ptheta in formula (51) [ibid]

   Map.phi1 <- vector("list",Dim)
   for(i in 1:Dim)
       { Map.phi1[[i]] <- function(x) phi@Map[[i]](x)* p(distr)(x)
         env.i <- environment(Map.phi1[[i]]) <- new.env()
         assign("i", i, envir=env.i)
         }

   phi1 <- EuclRandVariable(Map = Map.phi1, Domain = Reals())
##-t-## print(system.time({
   psi1 <- E(object=mu, fun = phi1, diagnostic=diagnostic)#, useApply = .useApply)
         if(diagnostic) diagn$psi1 <- attr(psi1,"diagnostic")
##-t-## }))

   ## obtaining IC psi  (formula (51))
   Map.psi <- vector("list",Dim)


##-t-##  print(system.time({
   for(i in 1:Dim)
     {

       env.i <- new.env()
       assign("i", i, envir=env.i)

       if(is(mu,"AbscontDistribution")){
            qxm <- q.l(mu)(x.mu.seq.b)

##-t-##  print(system.time({
            fct0.qq <- phi@Map[[i]](qxm)
##-t-##   }))
            fct0.q1 <-  rev(fct0.qq[iN.mu.1])
            fct0.q2 <-  rev(fct0.qq[iN.mu.2])
            fct0.q3 <-  rev(fct0.qq[iN.mu.3])
            phi0.q3 <-  h0.mu/100*rev(.csimpsum(fct0.q3))
            phi0.q2 <-  phi0.q3[1]+h0.mu*rev(.csimpsum(fct0.q2))
            phi0.q1 <-  phi0.q2[1]+h0.mu/100*rev(.csimpsum(fct0.q1))
            phi0.q <- c(phi0.q1,phi0.q2,phi0.q3)
            phi0a.q <- approxfun(x.mu.seq.a, phi0.q, yleft = 0, yright = rev(phi0.q)[1])
            psi0 <- function(x) {phi0a.q(p(mu)(x))-psi1[i]}
            assign("phi0.q", phi0.q, envir=env.i)
            assign("phi0a.q", phi0a.q, envir=env.i)
            assign("psi0", psi0, envir=env.i)
       }else{
            if(is(mu,"DiscreteDistribution")){
               phi.mu <- evalRandVar(phi,as.matrix(x.mu.seq))[i,,1]
               psi0.d <- cumsum(phi.mu*prob.mu)
               psi0.a <- approxfun(x.mu.seq, psi0.d, yleft = 0, yright = 0)
               psi0 <- function(x) psi0.a(x)
               assign("phi.mu", phi.mu, envir=env.i)
               assign("psi0.a", psi0.a, envir=env.i)
               assign("psi0.d", psi0.d, envir=env.i)
               assign("psi0", psi0, envir=env.i)
            }else{

               fct0 <- function(x,y) evalRandVar(phi, as.matrix(y))[i,,1]*(x<=y)
               diagn1 <- list()
               x.mu.seq.i <- quantile(seq(x.mu.seq),c(0,0.25,.5,.75,1))
               phi0 <- sapply(seq(x.mu.seq),
                           function(X){
                              fct <- function(y) fct0(x=x.mu.seq[X],y)
                              res <- do.call(E,c(list(object=mu, fun = fct,
                                                 useApply=.useApply),dotsInt))
                              if(diagnostic) if(X %in% x.mu.seq.i)
                                    diagn1[[paste(X)]] <<- attr(res,"diagnostic")
                              return(res)
                              })
               if(diagnostic) diagn[["phi0"]][[i]] <- diagn1

               phi0a <- approxfun(x.mu.seq, phi0, yleft = 0, yright = rev(phi0)[1])
               if(is(distr,"DiscreteDistribution"))
                     psi0 <- function(x) phi0a(x) * (x %in% support(mu))
               else  psi0 <- function(x) phi0a(x)

               assign("fct", fct, envir=env.i)
               assign("fct0", fct0, envir=env.i)
               assign("phi0", phi0, envir=env.i)
               assign("phi0a", phi0a, envir=env.i)
               assign("psi0", psi0, envir=env.i)
            }
       }

       Map.psi[[i]] <- psi0
       environment(Map.psi[[i]]) <- env.i

    }
##-t-##  }))
#   print(Map.psi)
   psi <-  EuclRandVariable(Map = Map.psi, Domain = Reals())

#   E2 <- E(object=distr, fun = psi %*%t(psi), useApply = .useApply)
   ## E2 = Cov_mu (psi)

   ### control: centering & standardization
   L2deriv.0 <- L2Fam@L2deriv[[1]]
##-t-##  print(system.time({
   E1 <- E(object=distr, fun = psi, diagnostic=diagnostic)
         if(diagnostic) diagn$E1 <- attr(E1,"diagnostic")

##-t-##  }))
##-t-##  print(system.time({
   E3 <- E(object=distr, fun = psi %*%t(L2deriv.0), diagnostic=diagnostic)
         if(diagnostic) diagn$E3 <- attr(E3,"diagnostic")
##-t-##  }))
   psi.0 <- psi - E1
   psi.01 <- as(distr::solve(E3)%*%psi.0,"EuclRandVariable")
   if(withplot)
      { for(i in 1:Dim)
         { dev.new()
           x0.mu.seq <- x.mu.seq
           if(is(mu,"AbscontDistribution")) x0.mu.seq <- q.l(mu)(x.mu.seq)
           plot(x0.mu.seq, sapply(x0.mu.seq,psi.01@Map[[i]]),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
         }}
##-t-##  print(system.time({
   E4 <- E(object=distr, fun = psi.01 %*%t(psi.01), diagnostic=diagnostic)
         if(diagnostic) diagn$E4 <- attr(E4,"diagnostic")
##-t-##  }))
   }
  E4 <- PosSemDefSymmMatrix(E4)

  psi <-  EuclRandVarList(psi.01)
  if(onedim){
     fctl <- vector("list",1)
     fct1 <- psi.01@Map[[1]]
     fctl[[1]] <- function(x) fct1(x)*liesInSupport(distr,x,checkFin=TRUE)
     env.1 <- environment(fctl[[1]]) <- new.env()
     assign("distr", distr, env.1)
     assign("fct1", fct1, env.1)
     assign("psi", psi, env.1)
     psi[[1]]@Map <- fctl
  }else{
     fctl <- vector("list",Dim)
     for(i in 1:Dim){
         fcti <-  psi[[1]]@Map[[i]]
         fctl[[i]] <- function(x) fcti(x)*liesInSupport(distr,x,checkFin=TRUE)
         env.i <- environment(fctl[[i]]) <- new.env()
         assign("i", i, env.i)
         assign("distr", distr, env.i)
         assign("fcti", fcti, env.i)
         assign("psi", psi, env.i)
     }
     psi[[1]]@Map <- fctl
  }
  nms <- names(c(main(param(L2Fam)),nuisance(param(L2Fam))))
  dimnames(E4) = list(nms,nms)
  if(diagnostic && !is.null(diagn)) class(diagn) <- "DiagnosticClass"
  if(diagnostic &&! withpreIC){
     attr(E4,"diagnostic") <- diagn
     if(!is.null(diagn)) class(attr(E4,"diagnostic"))<- "DiagnosticClass"
  }
  if(withpreIC) return(list(preIC=psi, Var=E4, diagnostic = diagn))
  else return(E4)
}


.oldCvMMDCovariance<- function(L2Fam, param, mu = distribution(L2Fam),
                            withplot = FALSE, withpreIC = FALSE,
                            N = getdistrOption("DefaultNrGridPoints")+1,
                            rel.tol=.Machine$double.eps^0.3,
                            TruncQuantile = getdistrOption("TruncQuantile"),
                            IQR.fac = 15,
                            ...){
   # preparations:

   N1 <- 2*N+1
   odd <- (1:N1)%%2==1

   param0 <- L2Fam@param
   dim0 <- dimension(param0)
#   print(param0)
   paramP <- param0
   paramP@main <- main(param)
   paramP@trafo <- diag(dim0)
#   print(paramP)
   L2Fam <- modifyModel(L2Fam, paramP)

#   print(L2deriv(L2Fam)[[1]]@Map)
   distr <- L2Fam@distribution

   ### get a sensible integration range:
   low0 <- q.l(distr)(TruncQuantile)
   up0 <- q.l(distr)(TruncQuantile, lower.tail = FALSE)
   m0 <- median(distr); s0 <- IQR(distr)
   low1 <- m0 - IQR.fac * s0
   up1  <- m0 + IQR.fac * s0
   low <- max(low0,low1); up <- min(up0,up1)

   ### get a sensible integration range:
   if(missing(mu)) mu <- distr
   low0.mu <- q.l(mu)(TruncQuantile)
   up0.mu <- q.l(mu)(TruncQuantile, lower.tail = FALSE)
   m0.mu <- median(mu); s0.mu <- IQR(mu)
   low1.mu <- m0.mu - IQR.fac * s0.mu
   up1.mu  <- m0.mu + IQR.fac * s0.mu
   low.mu <- max(low0.mu,low1.mu); up.mu <- min(up0.mu,up1.mu)

   if(is(distr,"DiscreteDistribution"))
       x.seq <-support(distr)
   else
       {if(is(distr,"AbscontDistribution")){
           x.seq0 <- seq(low, up, length = N1)
           h0 <- diff(x.seq0[1:2])
           x.seq <- x.seq0[odd]
          }else{
           x.seq <- seq(low,up, length = N)
          }
       }
   if(is(mu,"DiscreteDistribution"))
       x.mu.seq <- support(mu)
   else
       {if(is(mu,"AbscontDistribution")){
           x.mu.seq0 <- seq(low.mu, up.mu, length = N1)
           h0.mu <- diff(x.mu.seq0[1:2])
           x.mu.seq <- x.mu.seq0[odd]
          }else{
           x.mu.seq <- seq(low.mu, up.mu, length = N)
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
 #  print(Delta0)
   Delta1 <- approxfun(x.seq, Delta0, yleft = 0, yright = 0)
   if(is(distr,"DiscreteDistribution"))
      Delta <- function(x) Delta1(x) * (x %in% support(distr))
   else  Delta <- function(x) Delta1(x)
 #  print(Delta(x.seq))
 #  print(Delta(rnorm(100)))

   ## J = Var_Ptheta Delta
   J1 <- E(object=distr, fun = Delta)
#   print(J1)
   Delta.0 <- function(x) Delta(x) - J1
 #  print(Delta.0(x.seq))
 #  print(Delta.0(r(distr)(100))^2)
   #J <- distrExIntegrate(function(x) d(distr)(x)*Delta.0(x)^2, lower=low, upper=up)
   J <- E(object=distr, fun = function(x) Delta.0(x)^2 )
#   print(J)

   ### CvM-IC phi
   phi <- function(x) Delta.0(x)/J
#   print(head(sapply(x.seq,phi)))

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
 #  print(psi0)
   psi.1 <- approxfun(x.mu.seq, psi0, yleft = 0, yright = rev(psi0)[1])
   if(is(mu,"DiscreteDistribution"))
         psi <- function(x) (psi.1(x)-psi1) * (x %in% support(mu))
   else  psi <- function(x) psi.1(x)-psi1

#   print(sapply(x.mu.seq,psi))

   E2 <- E(object=distr, fun = function(x) psi(x)^2)
   L2deriv <- L2Fam@L2deriv[[1]]
   ## E2 = Cov_mu (psi)

#   ### control: centering & standardization
   E1 <- E(object=distr, fun = psi )
   E3 <- E(object=distr, fun = function(x) psi(x)*evalRandVar(L2deriv, x))
   psi.0 <- function(x) psi(x) - E1
   psi.01 <- function(x) psi.0(x)/E3
   if(withplot)
       { dev.new() #windows()
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
  # print("HLL")
  # print(x.seq0)
   for(i in 1:Dim)
       { if(is(distr,"AbscontDistribution")){
            #print(L2deriv@Map[[i]])
            fct <- NULL
            fct0 <- sapply(x.seq0, L2deriv@Map[[i]]) *
                           d(distr)(x.seq0)
            #print(fct0)
            Delta0 <-  h0*.csimpsum(fct0)
         }else{
            fct0 <- function(x,y) L2deriv@Map[[i]](x)*(x<=y)
            Delta0 <- sapply(x.seq, function(Y){ fct <- function(x) fct0(x,y=Y)
                                            return(E(object=distr, fun = fct))})
         }
         #print(Delta0)
         Delta1 <- approxfun(x.seq, Delta0, yleft = 0, yright = 0)
         if(is(distr,"DiscreteDistribution"))
               Delta <- function(x) Delta1(x) * (x %in% support(distr))
         else  Delta <- function(x) Delta1(x)
         Map.Delta[[i]] <- Delta
         env.i <- environment(Map.Delta[[i]]) <- new.env()
         assign("i", i, envir=env.i)
         assign("fct", fct, envir=env.i)
         assign("fct0", fct0, envir=env.i)
         assign("Delta", Delta, envir=env.i)
         assign("Delta0", Delta0, envir=env.i)
         assign("Delta1", Delta1, envir=env.i)
         if(withplot){
           dev.new()
           #windows()
           plot(x.seq, sapply(x.seq,Map.Delta[[i]]),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
         }

   }
   Delta <-  EuclRandVariable(Map = Map.Delta, Domain = Reals())



   ## J = Var_Ptheta Delta
##-t-##  print(system.time({
   J1 <- E(object=distr, fun = Delta)
##-t-##  }))
   Delta.0 <- Delta - J1
##-t-##  print(system.time({
   J <- E(object=distr, fun = Delta.0 %*%t(Delta.0))
##-t-##  }))
   ### CvM-IC phi
   phi <- as(distr::solve(J)%*%Delta.0,"EuclRandVariable")

   ## integrand phi x Ptheta in formula (51) [ibid]

   Map.phi1 <- vector("list",Dim)
   for(i in 1:Dim)
       { Map.phi1[[i]] <- function(x) evalRandVar(phi,x)[i] * p(distr)(x)
         env.i <- environment(Map.phi1[[i]]) <- new.env()
         assign("i", i, envir=env.i)
         }

   phi1 <- EuclRandVariable(Map = Map.phi1, Domain = Reals())
##-t-##  print(system.time({
   psi1 <- E(object=mu, fun = phi1)
##-t-##  }))


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
#       env.i0 <- environment(phi1) <- new.env()
#       assign("i", i, envir=env.i0)
       if(is(distr,"DiscreteDistribution"))
             psi0 <- function(x) phi0a(x) * (x %in% support(mu))
       else  psi0 <- function(x) phi0a(x)

       Map.psi[[i]] <- psi0
       env.i <- environment(Map.psi[[i]]) <- new.env()
       assign("i", i, envir=env.i)
       assign("fct", fct, envir=env.i)
       assign("fct01", fct0, envir=env.i)
       assign("fct0", fct0, envir=env.i)
       assign("psi0", psi0, envir=env.i)
       assign("phi0a", phi0a, envir=env.i)
       assign("phi0", phi0, envir=env.i)
    }
   psi <-  EuclRandVariable(Map = Map.psi, Domain = Reals())

##-t-##  print(system.time({
   E2 <- E(object=distr, fun = psi %*%t(psi))
##-t-##  }))
   ## E2 = Cov_mu (psi)

   ### control: centering & standardization
   L2deriv <- L2Fam@L2deriv[[1]]
##-t-##  print(system.time({
   E1 <- E(object=distr, fun = psi )
##-t-##  }))
##-t-##  print(system.time({
   E3 <- E(object=distr, fun = psi %*%t(L2deriv))
##-t-##  }))
   psi.0 <- psi - E1
   psi.01 <- as(distr::solve(E3)%*%psi.0,"EuclRandVariable")
   if(withplot)
      { for(i in 1:Dim)
         { devNew()
           plot(x.mu.seq, sapply(x.mu.seq,psi.01@Map[[i]]),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
         }}
##-t-##  print(system.time({
   E4 <- E(object=distr, fun = psi.01 %*%t(psi.01))
##-t-##  }))
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
P0 <- PoisFamily()
.oldCvMMDCovariance(P0,par=ParamFamParameter("lambda",1), withplot=TRUE)
.CvMMDCovariance(P0,par=ParamFamParameter("lambda",1), withplot=TRUE)
B0 <- BinomFamily(size=8, prob=0.3)
.oldCvMMDCovariance(B0,par=ParamFamParameter("",.3), withplot=TRUE)
.CvMMDCovariance(B0,par=ParamFamParameter("",.3), withplot=TRUE)
N0 <- NormLocationFamily();
.CvMMDCovariance(N0,par=ParamFamParameter("",0), withplot=TRUE, N = 200)
.oldCvMMDCovariance(N0,par=ParamFamParameter("",0), withplot=TRUE, N = 200)
C0 <- CauchyLocationFamily()
.oldCvMMDCovariance(C0,par=ParamFamParameter("",0), withplot=TRUE, N = 200)
.CvMMDCovariance(C0,par=ParamFamParameter("",0), withplot=TRUE, N = 200)
N1 <- NormScaleFamily()
.oldCvMMDCovariance(N1,par=ParamFamParameter("",1), withplot=TRUE, N = 200)
.CvMMDCovariance(N1,par=ParamFamParameter("",1), withplot=TRUE, N = 200)
NS <- NormLocationScaleFamily(); paramP <- ParamFamParameter(name = "locscale", main = c("loc"=0,"scale"=1),trafo = diag(2));
.oldCvMMDCovariance(NS,par=paramP, withplot=TRUE, N = 100)
.CvMMDCovariance(NS,par=paramP, withplot=TRUE, N = 100)
cls <- CauchyLocationScaleFamily();
.oldCvMMDCovariance(cls,par=ParamFamParameter("",0:1), withplot=TRUE, N = 200)
.CvMMDCovariance(cls,par=ParamFamParameter("",0:1), withplot=TRUE, N = 200)
Els <- L2LocationScaleFamily(loc = 0, scale = 1,
                  name = "Laplace Location and scale family",
                  centraldistribution = DExp(),
                  LogDeriv = function(x)  sign(x),
                  FisherInfo = diag(2),
                  trafo = diag(2))
.oldCvMMDCovariance(Els,par=ParamFamParameter("",0:1), withplot=TRUE, N = 100)
.CvMMDCovariance(Els,par=ParamFamParameter("",0:1), withplot=TRUE, N = 100)
Nb <- NbinomwithSizeFamily()
.oldCvMMDCovariance(Nb,par=ParamFamParameter(main=c(size=2.3,prob=0.3)), withplot=TRUE, N = 100)
.CvMMDCovariance(Nb,par=ParamFamParameter(main=c(size=2.3,prob=0.3)), withplot=TRUE, N = 100)

GF <- GammaFamily()
.oldCvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=4.3)), withplot=TRUE, N = 100)
.CvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=4.3)), withplot=TRUE, N = 100)
.oldCvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=0.3)), withplot=TRUE, N = 100)
.CvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=0.3)), withplot=TRUE, N = 100)

P0 <- PoisFamily()
B0 <- BinomFamily(size=8, prob=0.3)
N0 <- NormLocationFamily();
C0 <- CauchyLocationFamily()
cls <- CauchyLocationScaleFamily();
N1 <- NormScaleFamily()
NS <- NormLocationScaleFamily(); paramP <- ParamFamParameter(name = "locscale", main = c("loc"=0,"scale"=1),trafo = diag(2));
Els <- L2LocationScaleFamily(loc = 0, scale = 1,
                  name = "Laplace Location and scale family",
                  centraldistribution = DExp(),
                  LogDeriv = function(x)  sign(x),
                  FisherInfo = diag(2),
                  trafo = diag(2))
Nb <- NbinomwithSizeFamily()
GF <- GammaFamily()
system.time(print(.oldCvMMDCovariance(P0,par=ParamFamParameter("lambda",1))))
system.time(print(.CvMMDCovariance(P0,par=ParamFamParameter("lambda",1))))
system.time(print(.oldCvMMDCovariance(B0,par=ParamFamParameter("",.3))))
system.time(print(.CvMMDCovariance(B0,par=ParamFamParameter("",.3))))
system.time(print(.oldCvMMDCovariance(N0,par=ParamFamParameter("",0))))
system.time(print(.CvMMDCovariance(N0,par=ParamFamParameter("",0))))
system.time(print(.oldCvMMDCovariance(C0,par=ParamFamParameter("",0))))
system.time(print(.CvMMDCovariance(C0,par=ParamFamParameter("",0))))
system.time(print(.oldCvMMDCovariance(N1,par=ParamFamParameter("",1))))
system.time(print(.CvMMDCovariance(N1,par=ParamFamParameter("",1))))
system.time(print(.oldCvMMDCovariance(NS,par=paramP)))
system.time(print(.CvMMDCovariance(NS,par=paramP)))
system.time(print(.oldCvMMDCovariance(cls,par=ParamFamParameter("",0:1))))
system.time(print(.CvMMDCovariance(cls,par=ParamFamParameter("",0:1))))
system.time(print(.oldCvMMDCovariance(Els,par=ParamFamParameter("",0:1))))
system.time(print(.CvMMDCovariance(Els,par=ParamFamParameter("",0:1))))
system.time(print(.oldCvMMDCovariance(Nb,par=ParamFamParameter(main=c(size=2.3,prob=0.3)))))
system.time(print(.CvMMDCovariance(Nb,par=ParamFamParameter(main=c(size=2.3,prob=0.3)))))
system.time(print(.oldCvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=4.3)))))
system.time(print(.CvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=4.3)))))
system.time(print(.oldCvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=0.3)))))
system.time(print(.CvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=0.3)))))
system.time(print(.oldCvMMDCovariance(P0,par=ParamFamParameter("lambda",1),mu=Norm())))
system.time(print(.CvMMDCovariance(P0,par=ParamFamParameter("lambda",1),mu=Norm())))
system.time(print(.oldCvMMDCovariance(Nb,par=ParamFamParameter(main=c(size=2.3,prob=0.3)),mu=Norm())))
system.time(print(.CvMMDCovariance(Nb,par=ParamFamParameter(main=c(size=2.3,prob=0.3)),mu=Norm())))

" ## Code capsulated as string to avoid interpretation during R CMD Check

require(RobExtremes)
Pf <- ParetoFamily()
.oldCvMMDCovariance(Pf,par=param(Pf), withplot=TRUE)
.CvMMDCovariance(Pf,par=param(Pf), withplot=TRUE) ## better trust this one
GPDf <- GParetoFamily()
.oldCvMMDCovariance(GPDf,par=param(GPDf), withplot=TRUE)
.CvMMDCovariance(GPDf,par=param(GPDf), withplot=TRUE) ## better trust this one
GEVf <- GEVFamily()
.oldCvMMDCovariance(GEVf,par=param(GEVf), withplot=TRUE)
.CvMMDCovariance(GEVf,par=param(GEVf), withplot=TRUE) ## better trust this one
GEVuf <- GEVFamilyMuUnknown()
.oldCvMMDCovariance(GEVuf,par=param(GEVuf), withplot=TRUE)
.CvMMDCovariance(GEVuf,par=param(GEVuf), withplot=TRUE) ## better trust this one


system.time(print(.oldCvMMDCovariance(Pf,par=param(Pf))))
system.time(print(.CvMMDCovariance(Pf,par=param(Pf)))) ## better trust this one
system.time(print(.oldCvMMDCovariance(GPDf,par=param(GPDf))))
system.time(print(.CvMMDCovariance(GPDf,par=param(GPDf)))) ## better trust this one
system.time(print(.oldCvMMDCovariance(GEVf,par=param(GEVf))))
system.time(print(.CvMMDCovariance(GEVf,par=param(GEVf)))) ## better trust this one
system.time(print(.oldCvMMDCovariance(GEVuf,par=param(GEVuf))))
system.time(print(.CvMMDCovariance(GEVuf,par=param(GEVuf)))) ## better trust this one

"

set.seed(123)
x <- rnorm(100)
diF <- DiscreteDistribution(x)
system.time(print(.oldCvMMDCovariance(N0,par=ParamFamParameter("",0),mu=diF)))
system.time(print(.CvMMDCovariance(N0,par=ParamFamParameter("",0),mu=diF)))
system.time(print(.oldCvMMDCovariance(NS,par=paramP,mu=diF)))
system.time(print(.CvMMDCovariance(NS,par=paramP,mu=diF)))
IC <- .CvMMDCovariance(N0,par=ParamFamParameter("",0),withpreIC=TRUE)$preIC
sapply(c(-1e8,-3,0,0.1,1,2,1e8), function(x) evalRandVar(IC[[1]],x))
IC <- .CvMMDCovariance(GF,par=ParamFamParameter(main=c(scale=2.3,shape=4.3)),withpreIC=TRUE)$preIC
sapply(c(-1e8,-3,0,0.1,1,2,1e8), function(x) evalRandVar(IC[[1]],x))
IC <- .CvMMDCovariance(P0,par=ParamFamParameter("lambda",1),withpreIC=TRUE)$preIC
sapply(c(-1e8,-3,0,0.1,1,2,1e8), function(x) evalRandVar(IC[[1]],x))
IC <- .CvMMDCovariance(Nb,par=ParamFamParameter(main=c(size=2.3,prob=0.3)),withpreIC=TRUE)$preIC
sapply(c(-1e8,-3,0,0.1,1,2,1e8), function(x) evalRandVar(IC[[1]],x))
}
