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
A trace plot of the current evaluations of the approximate expected utility function. Separate lines are produced for the traces from Phases I and II of the ACE algorithm. These trace plots can be used to assess convergence. See Overstall & Woods (2017) for more details.
}
\references{
Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, in press. (\emph{arXiv:1501.00264}).

}
\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods & Maria Adamou
}
\note{For an example see \code{\link{ace}}.}

\seealso{
\code{\link{ace}}
}