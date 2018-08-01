################################################################
# return level - Plot functions in package distrMod
################################################################


setMethod("returnlevelplot", signature(x = "ANY",
                              y = "UnivariateDistribution"),
    function(x,    ### observations
             y,    ### distribution
             n = length(x), ### number of points to be plotted
             withIdLine = TRUE, ### shall line y=x be plotted in
             withConf = TRUE,   ### shall confidence lines be plotted
             withConf.pw  = withConf,   ### shall pointwise confidence lines be plotted
             withConf.sim = withConf,   ### shall simultaneous confidence lines be plotted
             plot.it = TRUE,    ### shall be plotted at all (inherited from stats::qqplot)
             datax = FALSE,     ### as in qqnorm
             MaxOrPOT = c("Max","POT"), ### used for block maxima or points over threshold
             npy = 365, ### number of observations per year
             threshold = if(is(y,"GPareto")) NA else 0,
             xlab = deparse(substitute(x)), ## x-label
             ylab = deparse(substitute(y)), ## y-label
             main = "",
             ...,                 ## further parameters
             width = 10,          ## width (in inches) of the graphics device opened
             height = 5.5,        ## height (in inches) of the graphics device opened}
             withSweave = getdistrOption("withSweave"), ## logical: if \code{TRUE}
             ##               (for working with \command{Sweave}) no extra device is opened and height/width are not set
             mfColRow = TRUE,     ## shall we use panel partition mfrow=c(1,1)?
             n.CI = n,            ## number of points to be used for CI
             with.lab = FALSE,     ## shall observation labels be plotted in
             lab.pts = NULL,      ## observation labels to be used
             which.lbs = NULL,    ## which observations shall be labelled
             which.Order = NULL,  ## which of the ordered (remaining) observations shall be labelled
             which.nonlbs = NULL, ## which of the non-labelled observations shall be plotted
             attr.pre = FALSE,    ## do indices refer to order pre or post ordering
             order.traf = NULL,   ## an optional trafo; by which the observations are ordered (as order(trafo(obs))
             col.IdL = "red",     ## color for the identity line
             lty.IdL = 2,         ## line type for the identity line
             lwd.IdL = 2,         ## line width for the identity line
             alpha.CI = .95,      ## confidence level
             exact.pCI = (n<100), ## shall pointwise CIs be determined with exact Binomial distribution?
             exact.sCI = (n<100), ## shall simultaneous CIs be determined with exact kolmogorov distribution?
             nosym.pCI = FALSE,   ## shall we use (shortest) asymmetric CIs?
             col.pCI = "orange",  ## color for the pointwise CI
             lty.pCI = 3,         ## line type for the pointwise CI
             lwd.pCI = 2,         ## line width for the pointwise CI
             pch.pCI = par("pch"),## symbol for points (for discrete mass points) in pointwise CI
             cex.pCI = par("cex"),## magnification factor for points (for discrete mass points) in pointwise CI
             col.sCI = "tomato2", ## color for the simultaneous CI
             lty.sCI = 4,         ## line type for the simultaneous CI
             lwd.sCI = 2,         ## line width for the simultaneous CI
             pch.sCI = par("pch"),## symbol for points (for discrete mass points) in simultaneous CI
             cex.sCI = par("cex"),## magnification factor for points (for discrete mass points) in simultaneous CI
             added.points.CI = TRUE, ## should the CIs be drawn through additional points?
             cex.pch = par("cex"),## magnification factor for the plotted symbols (for backward compatibility only, cex.pts in the sequel)
             col.pch = par("col"),## color for the plotted symbols (for backward compatibility only, col.pts in the sequel)
             cex.pts = 1,         ## magnification factor for labelled shown observations
             col.pts = par("col"),## color for labelled shown observations
             pch.pts = 19,        ## symbol for labelled shown observations
             cex.npts = 1,        ## magnification factor for non-labelled shown observations
             col.npts = grey(.5), ## color for non-labelled shown observations
             pch.npts = 20,       ## symbol for non-labelled shown observations
             cex.lbs = par("cex"),## magnification factor for the plotted observation labels
             col.lbs = par("col"),## color for the plotted observation labels
             adj.lbs = par("adj"),## adj parameter for the plotted observation labels
             alpha.trsp = NA,     ## alpha transparency to be added afterwards
             jit.fac = 0,         ## jittering factor used for discrete distributions
             jit.tol = .Machine$double.eps, ## tolerance for jittering: if distance 
                                 #is smaller than jit.tol, points are considered replicates
             check.NotInSupport = TRUE, ## shall we check if all x lie in support(y)
             col.NotInSupport = "red", ## if preceding check TRUE color of x if not in support(y)
             with.legend = TRUE,  ## shall a legend be plotted
             legend.bg = "white", ## background for the legend
             legend.pos = "topleft", ## position for the legend
             legend.cex = 0.8,     ## magnification factor for the legend
             legend.pref = "",     ## prefix for legend  text
             legend.postf = "",    ## postfix for legend text
             legend.alpha = alpha.CI, ## nominal level of CI
             debug = FALSE, ## shall additional debug output be printed out?
             withSubst = TRUE
    ){ ## return value as in stats::qqplot

    MaxOrPOT <- match.arg(MaxOrPOT)
    mc <- match.call(call = sys.call(sys.parent(1)))
    xcc <- as.character(deparse(mc$x))

   .mpresubs <- if(withSubst){
                   function(inx) 
                    .presubs(inx, c("%C", "%A", "%D" ),
                          c(as.character(class(x)[1]), 
                            as.character(date()), 
                            xcc))
               }else function(inx)inx

    if(missing(xlab)){mc$xlab <-  paste(gettext("Return level of"), xcc)}
    if(missing(ylab)){mc$ylab <-  gettext("Return period (years)")}
    if(missing(main)) mc$main <- gettext("Return level plot")
    mcl <- as.list(mc)[-1]
    mcl$datax <- NULL
    mcl$MaxOrPOT <- NULL
    mcl$npy <- NULL
    mcl$withSweave <- NULL
    mcl$mfColRow <- NULL
    mcl$type <-NULL
    mcl$debug <- NULL
    mcl$added.points.CI <- NULL
    if(is.null(mcl$datax)) datax <- FALSE
    force(x)

    thresh0 <- threshold 
    if(is(y,"GPareto")){ 
       if(is.na(threshold)) thresh0 <- location(y)
       y <- y - thresh0
       x <- x + thresh0
    }              

    rank0x <- rank(x)

    xj <- sort(x)

    if(any(.isReplicated(x, jit.tol))&&jit.fac>0)
       xj[.isReplicated(x, jit.tol)] <- jitter(x[.isReplicated(x, jit.tol)], factor=jit.fac)

    rank1x <- rank(xj)[rank0x]
    ind.x <- order(xj)
    xj <- sort(xj)

    p2rl <- function(pp){
               pp <- p(y)(pp)
               return(if(MaxOrPOT=="Max") -1/log(pp) else  1/(1-pp)/npy)
    }

    pp <- ppoints(length(xj))
    yc.o <- q.l(y)(pp)
    ycl <- p2rl(yc.o)

    ### extend range somewhat
#    pyn <- p(y)(10^(seq(-1, 3.75 + log10(npy), by = 0.1)))
    xyall <- force(sort(unique(c(yc.o,x,
                    q.l(y)(c(seq(0.01, 0.09, by = 0.01),(1:9)/10,
                         0.95, 0.99, 0.995, 0.999))
                         ))))
    rxyall  <- (max(xyall)-min(xyall))*0.6
    rxymean <- (max(xyall)+min(xyall))/2

    xyallc  <- seq(from=rxymean-rxyall,to=rxymean+rxyall, length.out=400)
#    print(xyallc)
    pxyall  <- p(y)(xyallc)
#    print(pxyall)

    pxyallc <- p2rl(xyallc)
     xyallc <-  xyallc[pxyall>0.00001 & pxyall<0.99999]
    pxyallc <- pxyallc[pxyall>0.00001 & pxyall<0.99999]

#    print(cbind(pxyallc,xyallc))

    if("support" %in% names(getSlots(class(y))))
       ycl <- sort(jitter(ycl, factor=jit.fac))

#-------------------------------------------------------------------------------
    alp.v <- .makeLenAndOrder(alpha.trsp,ind.x)
    alp.t <- function(x,a1) if(is.na(x)) x else addAlphTrsp2col(x,a1)
    alp.f <- if(length(alpha.trsp)==1L && is.na(alpha.trsp))
             function(x,a) x else function(x,a) mapply(x,alp.t,a1=a)

    if(missing(cex.lbs)) cex0.lbs <- par("cex")
    cex0.lbs <- .makeLenAndOrder(cex.lbs,ind.x)
    if(missing(adj.lbs)) adj0.lbs <- par("adj")
    adj0.lbs <- .makeLenAndOrder(adj.lbs,ind.x)
    if(missing(col.lbs)) col0.lbs <- par("col")
    col0.lbs <- alp.f(.makeLenAndOrder(col.lbs,ind.x),alp.v)
    if(missing(lab.pts)||is.null(lab.pts)) lab0.pts <- ind.x else
      lab0.pts <- .makeLenAndOrder(lab.pts,ind.x)

    lbprep <- .labelprep(x = x, y = yc.o[rank1x], lab.pts = lab0.pts,
                         col.lbs = col0.lbs, cex.lbs = cex0.lbs,
                         adj.lbs = adj0.lbs, which.lbs = which.lbs,
                         which.Order = which.Order, order.traf = order.traf,
                         which.nonlbs = which.nonlbs)

    n.ns <- length(lbprep$ns)
    n.s <- length(lbprep$ord)

    shown <- c(lbprep$ord,lbprep$ns)

    xs <- x[shown]
    ycs <- (ycl[rank1x])[shown]

    ordx <- order(xs)
    xso <- xs[ordx]
    ycso <- ycs[ordx]

    if(missing(cex.pch)) cex.pch <- par("cex")
    if(missing(col.pch)) col.pch <- par("col")
    if(missing(cex.pts)) cex.pts <- if(missing(cex.pch)) 1 else cex.pch
    if(missing(col.pts)) col.pts <- if(missing(col.pch)) par("col") else col.pch
    if(missing(pch.pts)) pch.pts <- 19
    if(missing(cex.npts)) cex.npts <- 1
    if(missing(col.npts)) col.npts <- par("col")
    if(missing(pch.npts)) pch.npts <- 20

    if(with.lab) lab.pts <- lbprep$lab.pts
    if(attr.pre){
       if(with.lab){
          col.lbs <- lbprep$col.lbs
          cex.lbs <- lbprep$cex.lbs
          adj.lbs <- lbprep$adj.lbs
       }
       cex.pts <- .makeLenAndOrder(cex.pts,ind.x)
       col.pts <- alp.f(.makeLenAndOrder(col.pts,ind.x),alp.v)
       pch.pts <- .makeLenAndOrder(pch.pts,ind.x)
       cex.pts <- cex.pts[shown]
       col.pts <- col.pts[shown]
       pch.pts <- pch.pts[shown]
    }else{
       ind.s <- 1:n.s
       ind.ns <- 1:n.ns
       if(with.lab){
          if(missing(lab.pts)||is.null(lab.pts)) lab.pts <- ind.ns else
             lab.pts <- .makeLenAndOrder(lab.pts,ind.ns)
          if(missing(cex.lbs)) cex.lbs <- par("cex")
          cex.lbs <- (.makeLenAndOrder(cex.lbs,ind.s))
          if(missing(adj.lbs)) adj.lbs <- par("adj")
          adj.lbs <- (.makeLenAndOrder(adj.lbs,ind.s))
          if(missing(col.lbs)) col.lbs <- par("col")
          col.lbs <- (alp.f(.makeLenAndOrder(col.lbs,ind.s),alp.v[lbprep$ord]))
       }
       cex.pts <- .makeLenAndOrder(cex.pts,ind.s)
       col.pts <- alp.f(.makeLenAndOrder(col.pts,ind.s),alp.v[lbprep$ord])
       pch.pts <- .makeLenAndOrder(pch.pts,ind.s)
       cex.npts <- .makeLenAndOrder(cex.npts,ind.ns)
       col.npts <- alp.f(.makeLenAndOrder(col.npts,ind.ns),alp.v[lbprep$ns])
       pch.npts <- .makeLenAndOrder(pch.npts,ind.ns)
       col.pts <- c(col.pts,col.npts)
       cex.pts <- c(cex.pts,cex.npts)
       pch.pts <- c(pch.pts,pch.npts)
    }
    cex.pts <- cex.pts[ordx]
    col.pts <- col.pts[ordx]
    pch.pts <- pch.pts[ordx]

#-------------------------------------------------------------------------------

    if(check.NotInSupport){
       xo <- xso #x[ord.x]
       nInSupp <- which(xo < q.l(y)(0))

       nInSupp <- unique(sort(c(nInSupp,which( xo > q.l(y)(1)))))
       if("support" %in% names(getSlots(class(y))))
          nInSupp <- unique(sort(c(nInSupp,which( ! xo %in% support(y)))))
       if("gaps" %in% names(getSlots(class(y))))
          nInSupp <- unique(sort(c(nInSupp,which( .inGaps(xo,gaps(y))))))
       if(length(nInSupp)){
#          col.pch[nInSupp] <- col.NotInSupport
          col.pts[nInSupp] <- col.NotInSupport
          if(with.lab)
#             col.lbs[ord.x[nInSupp]] <- col.NotInSupport
             col.lbs[nInSupp] <- col.NotInSupport
       }
    }

    if(n < length(x)){
       with.lab <- FALSE
       nos <- length(shown)
       idx <- sample(1:nos,size=n,replace=FALSE)
       cex.pts <- cex.pts[idx]
       col.pts <- col.pts[idx]
       pch.pts <- pch.pts[idx]
       xso <- xso[idx]
       ycso <- ycso[idx]
    }

    mcl <- .deleteItemsMCL(mcl)
    mcl$pch <- pch.pts
    mcl$cex <- cex.pts
    mcl$col <- col.pts

    mc$xlab <- .mpresubs(mcl$xlab)
    mc$ylab <- .mpresubs(mcl$ylab)

    if (!withSweave){
           devNew(width = width, height = height)
    }
    opar <- par("mfrow", no.readonly = TRUE)
    if(mfColRow) on.exit(do.call(par, list(mfrow=opar, no.readonly = TRUE)))

    if(mfColRow) opar1 <- par(mfrow = c(1,1), no.readonly = TRUE)

    ret <- list(x=xj,y=ycl)

    if(plot.it){
       xallc1 <- sort(c(xj,xyallc))
       yallc1 <- sort(c(ycl,pxyallc))
       mcl$x <- mcl$y <- NULL
       logs <- if(datax) "y" else "x"
       if(!is.null(mcl$log)){
           if(grepl("y", eval(mcl$log))) logs <- "xy"
           if(grepl("x",eval(mcl$log)))
              warning("The x axis is logarithmic anyway.")
           mcl$log <- NULL
       }
       if(datax){
          mcl$xlab <- mc$xlab
          mcl$ylab <- mc$ylab
          do.call(plot, c(list(x=xallc1, y=yallc1, log=logs,type="n"),mcl))
          do.call(points, c(list(x=xso, y=ycso), mcl))
       }else{
          mcl$ylab <- mc$xlab
          mcl$xlab <- mc$ylab
          do.call(plot,  c(list(x=yallc1, y=xallc1, log=logs,type="n"),mcl))
          do.call(points, c(list(x=ycso, y=xso), mcl))
       }
    }

    if(with.lab&& plot.it){
       lbprep$y0 <- p2rl(lbprep$y0)
       xlb0 <- if(datax) lbprep$x0 else lbprep$y0
       ylb0 <- if(datax) lbprep$y0 else lbprep$x0
       text(x = xlb0, y = ylb0, labels = lbprep$lab,
            cex = lbprep$cex, col = lbprep$col, adj = adj.lbs)
    }

    if(withIdLine){
       if(plot.it){
          if(datax){
             lines(xyallc,pxyallc,col=col.IdL,lty=lty.IdL,lwd=lwd.IdL)
          }else{
             lines(pxyallc,xyallc,col=col.IdL,lty=lty.IdL,lwd=lwd.IdL)
          }
       }
       qqb <- NULL
       if(#is(y,"AbscontDistribution")&&
       withConf){

          if(added.points.CI){
             xy <- unique(sort(c(x,xj,xyallc,yc.o)))
          }else{
             xy <- unique(sort(c(x,xj,yc.o)))
          }
          xy <- xy[!.NotInSupport(xy,y)]
          lxy <- length(xy)
          if(is(y,"DiscreteDistribution")){
             n0 <- min(n.CI, length(support(y)))
             n1 <- max(n0-lxy,0)
             if (n1 >0 ){
                 notyetInXY <- setdiff(support(y), xy)
                 xy0 <- sample(notyetInXY, n1)
                 xy <- sort(unique(c(xy,xy0)))
             }
          }else{
             if(lxy < n.CI){
                n1 <- (n.CI-lxy)%/%3
                xy0 <- seq(min(xy),max(xy),length=n1)
                xy1 <- r(y)(n.CI-lxy-n1)
                xy <- sort(unique(c(xy,xy0,xy1)))
             }
          }

        #qqb <- qqbounds(sort(unique(xy)),y,alpha.CI,n,withConf.pw, withConf.sim,
        #                   exact.sCI,exact.pCI,nosym.pCI, debug = debug)
        #qqb$crit <- p2rl(qqb$crit)
        if(plot.it){
          qqb <- .confqq(xy, y, datax, withConf.pw, withConf.sim, alpha.CI,
                      col.pCI, lty.pCI, lwd.pCI, pch.pCI, cex.pCI,
                      col.sCI, lty.sCI, lwd.sCI, pch.sCI, cex.sCI,
                  n, exact.sCI = exact.sCI, exact.pCI = exact.pCI,
                  nosym.pCI = nosym.pCI, with.legend = with.legend,
                  legend.bg = legend.bg, legend.pos = legend.pos,
                  legend.cex = legend.cex, legend.pref = legend.pref,
                  legend.postf = legend.postf, legend.alpha = legend.alpha,
                  qqb0=NULL, transf0=p2rl, debug = debug)
       }
    }}
    return(invisible(c(ret,qqb)))
    })

## into distrMod
setMethod("returnlevelplot", signature(x = "ANY",
                              y = "ProbFamily"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...){

    mc <- match.call(call = sys.call(sys.parent(1)))
    mcx <- as.character(deparse(mc$x))
    mcy <- as.character(deparse(mc$y))
    if(missing(xlab)) mc$xlab <- paste(gettext("Return Level of"), mcx)
    if(missing(ylab)) mc$ylab <- paste(gettext("Return Period at"), mcy)
    mcl <- as.list(mc)[-1]

    mcl$y <- yD <- y@distribution
    if(!is(yD,"UnivariateDistribution"))
       stop("Not yet implemented.")

    return(invisible(do.call(getMethod("returnlevelplot", signature(x="ANY", y="UnivariateDistribution")),
            args=mcl)))
    })

setMethod("returnlevelplot", signature(x = "ANY",
                              y = "Estimate"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...){

    mc <- match.call(call = sys.call(sys.parent(1)))
    mc1 <- match.call(call = sys.call(sys.parent(1)), expand.dots=FALSE)
    mcx <- as.character(deparse(mc$x))
    mcy <- as.character(deparse(mc$y))
    if(missing(xlab)) mc$xlab <- paste(gettext("Return Level of"), as.character(deparse(mc$x)))
    mcl <- as.list(mc)[-1]

    param <- ParamFamParameter(main=untransformed.estimate(y), nuisance=nuisance(y),
                               fixed=fixed(y))

    es.call <- y@estimate.call
    nm.call <- names(es.call)
    PFam <- NULL
    if("ParamFamily" %in% nm.call)
       PFam <- eval(as.list(es.call)[["ParamFamily"]])
    if(is.null(PFam))
       stop("There is no object of class 'ProbFamily' in the call of 'x'")

    PFam0 <- modifyModel(PFam, param)
    mcl$y <- PFam0
    if(missing(ylab)) mcl$ylab <- paste(gettext("Return Period at fitted"), name(PFam0), "\n -- fit by ", mcy)

    return(invisible(do.call(getMethod("returnlevelplot", signature(x="ANY", y="ProbFamily")),
            args=mcl)))
    })
