\name{plot.ace}
\alias{plot.ace}
\title{
Plot \code{ace} Objects
}
\description{
This function plots objects of class \code{"ace"}.
}
\usage{

\method{plot}{ace}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{"ace"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}
\value{
A trace plot of the current Monte Carlo evaluations to the expected utility function. Separate lines are produced for the traces from Phases I and II of the ACE algorithm. These trace plots can be used to assess convergence. See Overstall & Woods (2015) for more details.
}
\references{
Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.
}
\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk} & David C. Woods.
}
\note{For an example see \code{\link{ace}}.}

\seealso{
\code{\link{ace}}
}