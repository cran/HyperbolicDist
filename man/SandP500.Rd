\name{SandP500}
\alias{SandP500}
\title{Standard and Poor New York Stock Exchange market price index }
\description{
  This data set gives the value of the Standard and Poor New York
  Stock Exchange market price index at year end, from 1800 to 2001.
}
\usage{data(SandP500)}
\format{
  A vector of 202 observations.
}
\source{
  \url{http://www.globalfindata.com}
}
\references{
  Brown, Barry W., Spears, Floyd M. and Levy, Lawrence B. (2002)
  The log \emph{F}: a distribution for all seasons.
  \emph{Computational Statistics},
  \bold{17}, 47--58.
}
\examples{
data(SandP500)
## Consider proportional changes in the index
change<-SandP500[-length(SandP500)]/SandP500[-1]
hist(change)
## Fit hyperbolic distribution to changes
fit.hyperb(change)

}
\keyword{datasets}
