\name{calculate.range}
\alias{calculate.range}
\title{Range of a Hyperbolic Distribution}
\description{
  Given the parameter vector theta of a hyperbolic distribution,
  this function calculates the range over which the distribution
  is non-zero, to a specified tolerance. The parameterisation used
  is the pi/zeta one (see \code{\link{dhyperb}}). To use another
  parameterisation, use \code{\link{hyperb.change.pars}}. 
}
\usage{
calculate.range(theta, tol = 10^(-5))
}
\arguments{
  \item{theta}{Value of parameter vector specifying the hyperbolic
    distribution.}
  \item{tol}{Tolerance}
}
\details{
  The particular hyperbolic distribution being considered is specified
  by the value of the parameter value \code{theta}. The function calculates
  the effective range of the distribution, which is useful in calculating
  the distribution function and quantiles, but also in determining
  the range when plotting the density. By effective range is meant that
  the probability of an observation being greater than the upper end is
  less than the specified tolerance \code{tol}. Likewise for being smaller
  than the lower end of the range.
}
\value{
  A two-component vector giving the lower and upper ends of the range.
}
\references{ 
  Barndorff-Nielsen, O. and Blaesild, P (1983).
  Hyperbolic distributions. 
  In \emph{Encyclopedia of Statistical Sciences}, 
  eds., Johnson, N. L., Kotz, S. and Read, C. B., Vol. 3,
  pp. 700--707. New York: Wiley. 
}
\author{David Scott \email{d.scott@auckland.ac.nz},
  Jennifer Tso, Richard Trendall}

\seealso{
  \code{\link{dhyperb}}, \code{\link{hyperb.change.pars}}
}

\examples{
hyperb.range <- calculate.range(c(3,5,1,0),tol=10^(-2))
curve(dhyperb(x,c(3,5,1,0)),hyperb.range[1],hyperb.range[2])
}
\keyword{distribution}
