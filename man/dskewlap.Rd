\name{dskewlap}
\alias{dskewlap}
\title{Skew-Laplace Density}
\description{
  Calculates the density of the skew-Laplace distribution.
}
\usage{
dskewlap(x, theta)
}
\arguments{
  \item{x}{Vector of quantiles}
  \item{theta}{Vector of parameters of the skew-Laplace distribution:
    alpha, beta and mu}
}
\details{
  The central skew-Laplace has mode zero, and is a mixture of a (negative)
  exponential distribution with mean beta, and the negative of an 
  exponential distribution with mean alpha. The weights of the positive
  and negative components are proportional to their means.
  
  The general skew-Laplace distribution is a shifted central skew-Laplace
  distribution, where the mode is given by mu.
  
  The density is given by:
  \deqn{f(x)=\frac{1}{\alpha+\beta} e^{(x - \mu)/\alpha}}{ %
  f(x)=(1/(alpha+beta)) e^((x - mu)/alpha)}
  for \eqn{x\leq\mu}{x <= mu}, and
  \deqn{f(x)=\frac{1}{\alpha+\beta} e^{-(x - \mu)/\beta}}{ %
  f(x)=(1/(alpha+beta)) e^(-(x - mu)/beta)}
  for \eqn{x\geq\mu}{x >= mu}
}
\value{
  \code{dskewlap} gives the density of the skew-Laplace distribution.
}
\references{
  Fieller, N. J., Flenley, E. C. and Olbricht, W. (1992)
  Statistics of particle size data.
  \emph{Appl. Statist.},
  \bold{41}, 127--146.
}
\author{David Scott \email{d.scott@auckland.ac.nz}, 
Ai-Wei Lee, Richard Trendall}

\seealso{\code{\link{fit.hyperb}}}

\examples{
theta <- c(1,2,1)
curve(dskewlap(x,theta),from=-3,to=6,n=1000)
}
\keyword{distribution} 
