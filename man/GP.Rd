\name{GP}
\alias{distmat}
\alias{GPpred}
\alias{FisherScoring2D}

\title{
Functions for Fitting the Gaussian Process Emulator
}
\encoding{UTF-8}
\description{
These functions are used to fit a Gaussian process (GP) emulator for a stochastic function with a one-dimensional input space. The GP emulator is used within the approximate coordinate exchange (ACE) algorithm to approximate the expected utility function.
}
\usage{
distmat(Dij)

FisherScoring2D(par, Dij, z, dist, tol = 1e-10, max.iter = 25)

GPpred(paras, dist, z, newDij, Dij)
}

\arguments{
  \item{Dij}{
A vector of length \code{Q} specifying the values in the input space where the Monte Carlo integration approximation to the expected utility will be evaluated.
}
\item{paras}{
A vector of length 2 specifying the two parameters of the GP emulator (on the log scale). \code{paras[1]} specifies the correlation parameter and \code{paras[2]} specifies the nugget parameter.
}
\item{dist}{
A \code{Q} by \code{Q} matrix. The \eqn{ij}th element of \code{dist} specifies the squared difference of the ith and jth elements of \code{Dij}, i.e. a result of a call to \code{distmat}.
}
\item{z}{
A vector of length \code{Q} specifying the standardised evaluations of the Monte Carlo integration approximation to the expected utility at the elements of \code{Dij}.
}
\item{newDij}{
A vector specifying the values in the input space at which to approximate the expected utility.
}
\item{par}{
A vector of length 2 specifying the initial values of the GP emulator parameters (on the log scale) in the Fisher scoring algorithm.
}
\item{tol}{
Tolerance in the Fisher scoring algorithm. The default value is \code{1e-10}.
}
\item{max.iter}{
An integer specifying the maximum number of iterations of the Fisher scoring algorithm. The default value is \code{25}.
}
}
\details{
A key component in the ACE algorithm (see \code{\link{ace}}) is the GP emulator approximation to the expected utility function. This suite of functions fit the GP.

A GP emulator uses a correlation strcuture based on the distances, between the Monte Carlo evalutions of the expected utility function, in terms of the input space, to predict the value of the expected utility function.

The GP emulator has two unknown parameters: the so-called correlation and nugget parameters. These parameters are estimated via maximum likelihood using the Fisher scoring algorithm (Millar, 2011, pg 104). 

Once the parameters have been estimated, the GP emulator provides a predictive distribution of the expected utility function for any value in the input space. A point prediction is given by the predictive mean.

See Overstall & Woods (2015) for more details on the ACE algorithm and GP emulators.
}

\value{
\code{distmat} will return a \code{Q} by \code{Q} matrix where the \eqn{ij}th element is the squared difference between the \eqn{i}th and \eqn{j}th elements of the argument \code{Dij}. Note that the elements on the leading diagonal will be zero.

\code{FisherScoring2D} will return a list with two elements: \code{par} - the maximum likelihood estimates (on the log scale) of the two parameters of the GP emulator; and \code{singular} - this returns \code{1} if the final maximum likelihood estimates result in the GP correlation matrix being computationally singular, and \code{0} otherwise.

\code{GPpred} will return a vector of the same length as the argument \code{newDij} whose elements are the predictive means of the standardised output of the expected utility function.
}

\references{
Millar, R. (2011). \emph{Maximum Likelihood Estimation and Inference}. Wiley.

Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.
}

\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk} & David C. Woods.
}

\note{
Note that a typical user will not call these functions directly.
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

Dij<-sort(as.vector(randomLHS(n=20,k=1)))*2-1
## Generate a vector at which to evaluate the noisy evaluation of the expected utilty

neval<-noisyexpectutil(utility=utility,B=1000,d=d,i=3,j=4,Dij=Dij)
## Evaluate the noisy evaluation of the expected utility at each value of Dij, where the 
## 3rd row and 4th column is being varied.

neval

#$z
# [1] -17.93285 -18.23384 -18.03137 -17.95596 -17.98777 -18.16553 -18.01260
# [8] -18.12832 -18.40023 -18.22158 -18.25258 -18.34298 -18.24723 -18.58172
#[15] -18.74440 -18.67967 -18.61615 -19.00058 -18.69601 -19.11678
#
#$Dij
# [1] -0.92429129 -0.86746478 -0.72226793 -0.62888788 -0.51356605 -0.46609271
# [7] -0.39157531 -0.21246787 -0.16662251 -0.05653405  0.03999944  0.18921983
#[13]  0.21216919  0.33466835  0.43899895  0.54763512  0.67125147  0.78394404
#[19]  0.89606180  0.92026923

## EXAMPLE OF distmat

dist<-distmat(Dij)

dist[1:4,1:4]
## Print out first 4 rows and columns.

#            [,1]        [,2]        [,3]        [,4]
#[1,] 0.000000000 0.003229251 0.040813436 0.087263173
#[2,] 0.003229251 0.000000000 0.021082127 0.056918941
#[3,] 0.040813436 0.021082127 0.000000000 0.008719834
#[4,] 0.087263173 0.056918941 0.008719834 0.000000000

## EXAMPLE OF FisherScoring2D

mz<-mean(neval$z)
sz<-sd(neval$z)
z<-(neval$z - mz)/sz
## Standardise neval$z and call this z

FS<-FisherScoring2D(par = c(0,0), Dij = Dij, z = z, dist = dist)
## Find MLEs of GP parameters (on log scale). 

FS

#$par
#[1] -1.8357336 -0.4662313

## EXAMPLE OF GPpred

newDij<-sort(as.vector(randomLHS(n=5,k=1)))*2-1
newneval<-noisyexpectutil(utility=utility,B=1000,d=d,i=3,j=4,Dij=newDij)
## Generate a vector of points to predict the expected utility.

mz + sz*GPpred(para = FS$par, dist = dist, z = z, newDij = newDij, Dij = Dij)
## Predictions from GP emulator.
#[1] -18.04710 -18.09884 -18.27646 -18.42334 -18.71939

newneval$z
## Values generated using Monte Carlo integration.
#[1] -17.80276 -17.90503 -18.36726 -18.40718 -18.78144
}
