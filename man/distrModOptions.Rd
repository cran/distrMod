\name{distrModOptions}
\alias{distrModOptions}
\alias{getdistrModOption}
\alias{distrModoptions}
\alias{show.details}

\title{Function to change the global variables of the package `distrMod' }
\description{With \code{distrModOptions} you can inspect and change 
the global variables of the package \pkg{distrMod}. }
\usage{
distrModOptions(...)
getdistrModOption(x)
distrModoptions(...)
}
\arguments{
  \item{\dots}{ any options can be defined, using name = value or by passing a 
  list of such tagged values. }
  \item{x}{ a character string holding an option name.}
}
\details{Invoking \code{distrModoptions()} with no arguments returns a list 
 with the current values of the options.  To access the value of a single 
 option, one should use \code{getdistrModOption("show.details")}, 
 e.g., rather than \code{distrModoptions("show.details")} which is 
 a \emph{list} of length one.}
\value{
\code{distrModoptions()} returns a list of the global options of 
\pkg{distrMod}. \cr
\code{distrModoptions("show.details")} returns  the global option 
\code{show.details} as a list of length 1. \cr
\code{distrModoptions("show.details" = "minimal")} sets the value of the 
global option  \code{show.details} to \code{"minimal"}. 
\code{getdistrModOption("show.details")} the current value set for 
option \code{show.details}. 
}
\section{distrModoptions}{
For compatibility with spelling in package \pkg{distr}, \code{distrModoptions} is 
just a synonym to \code{distrModoptions}.
}

\section{Currently available options}{
\describe{
  \item{show.details}{degree of detailedness for method \code{show} for
  objects of classes of the \pkg{distrXXX} family of packages. Possible
  values are
  \describe{
  \item{\code{"maximal"}}{all information is shown}
  \item{\code{"minimal"}}{only the most important information is shown}
  \item{\code{"medium"}}{somewhere in the middle; see actual \code{show}-methods
  for details.}
  }
  The default value is \code{"maximal"}.}
}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
        Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[base]{options}}, \code{\link[base:options]{getOption}},
\code{\link[distr]{distroptions}}, \code{\link[distr:distroptions]{getdistrOption}}}
\examples{
distrModoptions()
distrModoptions("show.details")
distrModoptions("show.details" = "maximal")
distrModOptions("show.details" = "minimal")
# or
getdistrModOption("show.details")
}
\keyword{misc}
\keyword{distribution}
\concept{global options}
\concept{options}
