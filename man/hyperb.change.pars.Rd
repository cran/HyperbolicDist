\name{hyperb.change.pars}
\alias{hyperb.change.pars}

\title{Change Parameterisations of the Hyperbolic Distribution}

\description{
	This function interchanges between the following 3 parameterisations
	of the Hyperbolic Distribution:
  1. pi, zeta, delta, mu
  2. alpha, beta, delta, mu
  3. phi, gamma, delta, mu}

\usage{
hyperb.change.pars(from, to, theta, no.names=FALSE)
}

\arguments{
  \item{from}{the set of parameters to change from}
  \item{to}{the set of parameters to change to}
  \item{theta}{"from" parameter vector consisting of 4 numerical elements}
  \item{no.names}{logical; when \code{TRUE}, suppresses the parameter
		\code{names} in the output}
}

\details{
In the 3 parameterisations, the following must be positive:
  1. zeta, delta
  2. alpha, delta
  3. phi, gamma, delta
	Furthermore, note that in the second parameterisation alpha must must
	be greater than the absolute value of beta.
}

\value{
	A numerical vector of length 4 representing \code{theta} in the
	\code{to} parameterisation.
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

\seealso{\code{\link{dhyperb}}}

\examples{
theta1 <- c(-2,1,3,0) # Parameter vector in parameterisation 1
theta2 <- hyperb.change.pars(1,2,theta1) # Convert to parameterisation 2
theta2 # Parameter vector in parameterisation 2
hyperb.change.pars(2,1,as.numeric(theta2)) # Convert back to parameterisation 1
}

\keyword{distribution}
