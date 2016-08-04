\name{pval}
\alias{pval}

\title{
Bayesian Hypothesis Test
}
\encoding{UTF-8}
\description{
This function performs a Bayesian hypothesis test to determine the posterior probability that the proposed design has a larger expected utility than the current design.}
\usage{
pval(oldeval, neweval, binary)
}

\arguments{
  \item{oldeval}{
A vector of length B containing the value of the utility function, under the current design, evaluated at B values generated from the joint distribution of model parameters and unobserved responses.
}
  \item{neweval}{
A vector of length B containing the value of the utility function, under the proposed design, evaluated at B values generated from the joint distribution of model parameters and unobserved responses.
}
  \item{binary}{
A logical argument indicating whether the utility function outputs binary or real numbers. In some cases, the utility function is an indicator function of some event. The expected utility function will then be the expected posterior probability of the event. Utility functions such as Shannon information gain and negative squared error loss output real numbers.
}
}
\details{
In the course of the ACE algorithm, we will need to make pairwise comparisons of designs, where we compare according to the expectation of the utility function. We cannot evaluate the expectation of the utility function but we can generate a sample from a distribution whose theoretical mean is equal to the expected utility.

We can compare two designs (current and proposed) by generating two samples from the respective distributions. The posterior probability that the proposed design has a larger value of the expected utility can be given by the a Bayesian hypothesis test.

See Overstall & Woods (2015) for more details.
}

\value{
The function will return a scalar giving the posterior probability that the proposed design has a larger value of the expected utility.
}

\references{
Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.
}

\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk} & David C. Woods.
}

\note{
Note that a typical user will not call this function directly.
}

\examples{
set.seed(1)
## Set seed for reproducibility

oldx <- cbind( 1 , runif(10) )
newx <- cbind( 1 , runif(10) )
## Create model matrices for two designs.

beta <- cbind( runif(10) , runif(10) )
oldeval<-utilglm(x = oldx, beta = beta, family = binomial(link = "probit"))
beta <- cbind( runif(10) , runif(10) )
neweval<-utilglm(x = newx, beta = beta, family = binomial(link = "probit"))

## Evaluate utilty function corresponding to Bayesian D-optimality
## for each of the 10 values in the prior sample for each design.

pval(oldeval = oldeval, neweval = neweval, binary = FALSE)
## Calculate the posterior probability from the Bayesian hypothesis test:
# 0.005829958 
}
