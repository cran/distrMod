\name{qqplot}
\docType{methods}
\title{Methods for Function qqplot in Package `distrMod'}
\usage{
qqplot(x, y, ...)
\S4method{qqplot}{ANY,UnivariateDistribution}(x,y,
    n = length(x), withIdLine = TRUE,
    withConf = TRUE, withConf.pw  = withConf, withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),ylab = deparse(substitute(y)),
    ..., width = 10, height = 5.5, withSweave = getdistrOption("withSweave"),
    mfColRow = TRUE, n.CI = n, withLab = FALSE, lab.pts = NULL, which.lbs = NULL,
    which.Order = NULL, order.traf = NULL,
    col.IdL = "red", lty.IdL = 2, lwd.IdL = 2, alpha.CI = .95,
    exact.pCI = (n<100), exact.sCI = (n<100), nosym.pCI = FALSE,
    col.pCI = "orange", lty.pCI = 3, lwd.pCI = 2, pch.pCI = par("pch"),
    cex.pCI = par("cex"),
    col.sCI = "tomato2", lty.sCI = 4, lwd.sCI = 2, pch.sCI = par("pch"),
    cex.sCI = par("cex"),
    cex.pch = par("cex"), col.pch = par("col"),
    cex.lbl = par("cex"), col.lbl = par("col"), adj.lbl = NULL,
    jit.fac = 0, check.NotInSupport = TRUE, col.NotInSupport = "red",
    with.legend = TRUE, legend.bg = "white",
    legend.pos = "topleft", legend.cex = 0.8)
\S4method{qqplot}{ANY,ProbFamily}(x, y,
   n = length(x), withIdLine = TRUE, withConf = TRUE,
   withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...)
}
\alias{qqplot}
\alias{qqplot-methods}
\alias{qqplot,ANY,ProbFamily-method}
\alias{qqplot,ANY,UnivariateDistribution-method}

\arguments{
\item{x}{data to be checked for compatibility with distribution/model \code{y}.}
\item{y}{object of class \code{"UnivariateDistribution"} or of
class \code{"ProbFamily"}.}
\item{n}{numeric; assumed sample size (by default length of \code{x}).}
\item{withIdLine}{logical; shall line \code{y = x} be plotted in?}
\item{withConf}{logical; shall confidence lines be plotted?}
\item{withConf.pw}{logical; shall pointwise confidence lines be plotted?}
\item{withConf.sim}{logical; shall simultaneous confidence lines be plotted?}
\item{plot.it}{logical; shall be plotted at all (inherited from
\code{\link[stats:qqnorm]{qqplot}})?}
\item{xlab}{x-label}
\item{ylab}{y-label}
\item{\dots}{further parameters for method \code{qqplot} with signature
\code{ANY,UnivariateDistribution} or with function \code{plot}}
\item{width}{width (in inches) of the graphics device opened}
\item{height}{height (in inches) of the graphics device opened}
\item{withSweave}{logical: if \code{TRUE} (for working with \command{Sweave})
        no extra device is opened and height/width are not set}
\item{mfColRow}{shall default partition in panels be used --- defaults to \code{TRUE}}
\item{n.CI}{numeric; number of points to be used for confidence interval}
\item{withLab}{logical; shall observation labels be plotted in?}
\item{lab.pts}{character or \code{NULL}; observation labels to be used}
\item{which.lbs}{integer or \code{NULL}; which observations shall be labelled}
\item{which.Order}{integer or \code{NULL}; which of the ordered (remaining)
observations shall be labelled}
\item{order.traf}{function or \code{NULL}; an optional trafo by which the
observations are ordered (as order(trafo(obs)).}
\item{col.IdL}{color for the identity line}
\item{lty.IdL}{line type for the identity line}
\item{lwd.IdL}{line width for the identity line}
\item{alpha.CI}{confidence level}
\item{exact.pCI}{logical; shall pointwise CIs be determined with exact
Binomial distribution?}
\item{exact.sCI}{logical; shall simultaneous CIs be determined with
exact Kolmogorov distribution?}
\item{nosym.pCI}{logical; shall we use (shortest) asymmetric CIs?}
\item{col.pCI}{color for the pointwise CI}
\item{lty.pCI}{line type for the pointwise CI}
\item{lwd.pCI}{line width for the pointwise CI}
\item{pch.pCI}{symbol for points (for discrete mass points) in pointwise CI}
\item{cex.pCI}{magnification factor for points (for discrete mass points) in
pointwise CI}
\item{col.sCI}{color for the simultaneous CI}
\item{lty.sCI}{line type for the simultaneous CI}
\item{lwd.sCI}{line width for the simultaneous CI}
\item{pch.sCI}{symbol for points (for discrete mass points) in simultaneous CI}
\item{cex.sCI}{magnification factor for points (for discrete mass points) in
simultaneous CI}
\item{cex.pch}{magnification factor for the plotted symbols}
\item{col.pch}{color for the plotted symbols}
\item{cex.lbl}{magnification factor for the plotted observation labels}
\item{col.lbl}{color for the plotted observation labels}
\item{adj.lbl}{adj parameter for the plotted observation labels}
\item{jit.fac}{jittering factor used for discrete distributions}
\item{check.NotInSupport}{logical; shall we check if all \code{x}-quantiles lie
in \code{support(y)}?}
\item{col.NotInSupport}{logical; if preceding check \code{TRUE} color of
\code{x}-quantiles if not in \code{support(y)}}
\item{with.legend}{logical; shall a legend be plotted?}
\item{legend.bg}{background color for the legend}
\item{legend.pos}{position for the legend}
\item{legend.cex}{magnification factor for the legend}
}
\description{
  We generalize function \code{\link[stats:qqnorm]{qqplot}} from package \pkg{stats} to
  be applicable to distribution and probability model objects. In this context,
  \code{qqplot} produces a QQ plot of data (argument \code{x}) against
   a (model) distribution.
   Graphical parameters may be given as arguments to \code{qqplot}.
}
\value{
    As for function \code{\link[stats:qqnorm]{qqplot}} from package \pkg{stats}: a
  list with components
\item{x}{The x coordinates of the points that were/would be plotted}
  \item{y}{The corresponding quantiles of the second distribution,
           \emph{including \code{\link{NA}}s}.}
}
\details{
\describe{
\item{qqplot}{\code{signature(x = "ANY", y = "UnivariateDistribution")}:
produces a QQ plot of a dataset \code{x} against the theoretical
quantiles of distribution \code{y}.}
\item{qqplot}{\code{signature(x = "ANY", y = "ProbFamily")}:
produces a QQ plot of a dataset \code{x} against the theoretical
quantiles of the model distribution of model \code{y}. Passed through
the \code{\dots} argument, all arguments valid for
\code{signature(x = "ANY", y = "UnivariateDistribution")}
are also valid for this signature.}
}
}

\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language}.
  Wadsworth & Brooks/Cole.
}
\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
}
\seealso{
  \code{\link[stats:qqnorm]{qqplot}} from package \pkg{stats} -- the standard QQ plot
  function,  \code{\link[distr]{qqplot}} from package \pkg{distr} for
  comparisons of distributions, and
  \code{\link{qqbounds}}, used by \code{qqplot} to produce confidence
   intervals.
}
\examples{
qqplot(r(Norm(15,sqrt(30)))(40), Chisq(df=15))
qqplot(r(Norm(15,sqrt(30)))(40), NormLocationFamily())
}
\keyword{hplot}
\keyword{distribution}
