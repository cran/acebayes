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

If the object is a result of a direct call to \code{\link{aceglm}} or \code{\link{acenlm}}, then the argument \code{criterion} will be printed, otherwise the statement \code{User-defined utility} will be printed.

Also printed are the number of runs, factors, Phase I and II iterations of the ACE algorithm and the computational time required.

For more details on the ACE algorithm, see Overstall & Woods (2017).
}

\references{
Overstall, A.M. & Woods, D.C. (2016). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, in press. (\emph{arXiv:1501.00264}).
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods & Maria Adamou
}

\note{For examples see \code{\link{ace}}, \code{\link{aceglm}}, and \code{\link{acenlm}}.}

\seealso{
\code{\link{ace}}, \code{\link{aceglm}}, \code{\link{acenlm}}.
}
