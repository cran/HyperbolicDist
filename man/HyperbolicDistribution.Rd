\name{HyperbolicDistribution}
\alias{hyperb}

\title{
  The Library `HyperbolicDist': Summary Information
}

\description{
This library provides a collection of functions for working with
the hyperbolic distribution. Included are the 
density function, distribution function, quantiles, 
random number generation and fitting the hyperbolic distribution
to data (\code{fit.hyperb}). The function \code{hyperb.change.pars}
will interchange parameter values between different parameterisations.
The mean, variance and mode of a given hyperbolic distribution are
given by \code{hyperb.mean}, \code{hyperb.var}, and \code{hyperb.mode}
respectively. For assessing the fit of the hyperbolic distribution to
a set of data, the log-histogram is useful. See \code{\link{log.hist}}.
}

\seealso{
\code{\link{dhyperb}}, \code{\link{fit.hyperb}}, 
\code{\link{log.hist}}, \code{\link{hyperb.change.pars}}, and
\code{\link{hyperb.mean}}.
}

\author{David Scott \email{d.scott@auckland.ac.nz}}

\section{Acknowledgements}{
  A number of students have worked on the package: Ai-Wei Lee, Jennifer Tso,
  and Richard Trendall.
  
  Thanks to Ross Ihaka and Paul Murrell for their willingness to answer
  my questions, and to all the core group for the development of R.
}
\references{
  Barndorff-Nielsen, O. (1977)
  Exponentially decreasing distributions for the logarithm of particle size,
  \emph{Proc. Roy. Soc. Lond.},
  A\bold{353}, 401--419.
  
  Barndorff-Nielsen, O. and Blaesild, P (1983).
  Hyperbolic distributions. 
  In \emph{Encyclopedia of Statistical Sciences}, 
  eds., Johnson, N. L., Kotz, S. and Read, C. B., Vol. 3,
  pp. 700--707. New York: Wiley.
  
  Fieller, N. J., Flenley, E. C. and Olbricht, W. (1992)
  Statistics of particle size data.
  \emph{Appl. Statist.},
  \bold{41}, 127--146.
}


\section{LICENCE}{
This library and its documentation are usable under the terms of the 
"GNU General Public License", a copy of which is distributed with the 
package.
}


\keyword{distribution}
