\name{hyperb.mean}
\alias{hyperb.mean}
\alias{hyperb.var}
\alias{hyperb.mode}
\alias{R.lambda}
\alias{S.lambda}
\title{Mean, Variance and Mode of the Hyperbolic Distribution}
\description{
  Functions to calculate the mean, variance, and mode of a specific
  hyperbolic distribution. \code{R.lambda} and \code{S.lambda} are
  functions derived from Bessel functions which are used in the
  calculation of the mean and variance.
}
\usage{
hyperb.mean(theta)
hyperb.var(theta)
hyperb.mode(theta)
R.lambda(zeta, lambda = 1)
S.lambda(zeta, lambda = 1)

}
\arguments{
  \item{theta}{Parameter vector of the hyperbolic distribution}
  \item{zeta}{Value of the parameter zeta of the hyperbolic distribution}
  \item{lambda}{Parameter related to order of Bessel functions}
}

\value{
  \code{hyperb.mean} gives the mean of the hyperbolic distribution,
  \code{hyperb.var} the variance, and \code{hyperb.mode} the mode. 
  The formulae used are as given in Barndorff-Nielsen and Blaesild 
  (1983), p. 702. The functions \code{R.lambda} and \code{S.lambda}
  are used in the calculation of the mean and variance. They are 
  functions of the Bessel functions of the third kind,
  implemented in R as \code{\link{besselK}}.
  
  The parameterisation of the hyperbolic distribution used for this
  and other components of the \code{hyperbolic} package is the
  pi/zeta one. See \code{hyperb.change.pars} to transfer between
  parameterisations.
}
\references{
  Barndorff-Nielsen, O. and Blaesild, P (1983).
  Hyperbolic distributions. 
  In \emph{Encyclopedia of Statistical Sciences}, 
  eds., Johnson, N. L., Kotz, S. and Read, C. B., Vol. 3,
  pp. 700--707. New York: Wiley.
}

\author{David Scott \email{d.scott@auckland.ac.nz}, Richard Trendall}

\seealso{
  \code{\link{dhyperb}}, \code{\link{hyperb.change.pars}}, 
  \code{\link{besselK}}
}

\examples{
theta <- c(2,2,2,2)
hyperb.mean(theta)
hyperb.var(theta)
hyperb.mode(theta)
}
\keyword{distribution}
