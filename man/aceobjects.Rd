\name{aceobjects}
\alias{print.ace}
\alias{summary.ace}

\title{
Print and Summary of \code{ace} Objects
}
\description{
These function prints and summarises objects of class \code{"ace"}.
}
\usage{
\method{print}{ace}(x, ...)
\method{summary}{ace}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"ace"}.
}
  \item{object}{
An object of class \code{"ace"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
These functions both provide a print out with the following information.

If the object is a result of a direct call to \code{\link{aceglm}}, then the pesudo-Bayesian criterion will be printed, otherwise the statment \code{User-defined utility} will be printed.

Also printed are the number of runs, factors, Phase I and II iterations of the ACE algorithm and the computational time required.

For more details on the ACE algorithm, see Overstall & Woods (2015).
}

\references{
Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.
}

\author{
Antony M. Overstall (\email{Antony.Overstall@glasgow.ac.uk}) & David C. Woods.
}

\note{For examples see \code{\link{ace}} and \code{\link{aceglm}}.}

\seealso{
\code{\link{ace}}, \code{\link{aceglm}}.
}
