\name{noisyexpectutil}
\alias{noisyexpectutil}

\title{
Noisy Evaluation of Expected Utility Function
}
\encoding{UTF-8}
\description{
This function calculates a noisy evaluation of the expected utility function, using Monte Carlo integration, for a specified element (coordinate) of the current design.
}
\usage{
noisyexpectutil(utility,B,d,i,j,Dij)
}

\arguments{
  \item{utility}{
A function with two arguments: \code{d} and \code{B}. It should return a vector of length \code{B} where each element gives the value of the utility function for design \code{d}, for a value generated from the joint distribution of model parameters and unobserved responses.
}
\item{B}{
An integer giving the size of the sample generated from the joint distribution of model parameters and unobserved responses.
}
\item{d}{
An \code{n} by \code{k} matrix giving the current design.
}
\item{i}{
The row number of \code{d} identifying the element of \code{d} to vary. Should be in 1,...,n. 
}
\item{j}{
The column number of \code{d} identifying the element of \code{d} to vary. Should be in 1,...,k. 
}
\item{Dij}{
A vector of length Q giving values at which to evaluate the noisy evaluation of the expected utility. 
}
}
\details{
A key component in the ACE algorithm (see \code{\link{ace}}) is the GP emulator approximation to the expected utility function.

A vector of length \code{Q} of Monte Carlo integration evaluations of the expected utility function are required to fit the GP emulator.

Given a utility function as an argument, this function will return the required vector.

For more details on the ACE algorithm, see Overstall & Woods (2015).
}

\value{
The function will return a list with two elements.
  \item{Dij}{
The argument \code{Dij}. 
}
  \item{z}{
A vector of length Q giving the noisy evaluations of the expected utility at each value of \code{Dij}. 
}
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

\seealso{
\code{\link[=distmat]{GP}}
}
\examples{
utility<-function(d, B){
x<-cbind(1,d)
beta<-cbind(runif(n=B,min=-3,max=3),
runif(n=B,min=4,max=10),
runif(n=B,min=5,max=11),
runif(n=B,min=-6,max=0),
runif(n=B,min=-2.5,max=3.5))
utilglm(x=x, beta=beta, family=binomial, criterion="D")}
## Create a utility function for logistic regression under Bayesian D-optimality 

set.seed(1)
## Set seed for reproducibility

d<-randomLHS(n=6,k=4)*2-1
## Generate a current design

Dij<-sort(as.vector(randomLHS(n=5,k=1)))*2-1
## Generate a vector at which to evaluate the noisy evaluation of the expected utilty

neval<-noisyexpectutil(utility=utility,B=1000,d=d,i=3,j=4,Dij=Dij)
## Evaluate the noisy evaluation of the expected utility at each value of Dij, where the 
## 3rd row and 4th column is being varied.

neval

#$z
#[1] -17.89424 -18.18582 -18.15289 -18.49041 -18.64084
#
#$Dij
#[1] -0.8734913 -0.5602135 -0.1717284  0.4074537  0.6979189
}
