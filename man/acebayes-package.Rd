\name{acebayes-package}
\alias{acebayes-package}
\alias{acebayes}
\docType{package}
\title{
Optimal Bayesian Experimental Design using the Approximate Coordinate Exchange (ACE) Algorithm
}
\description{
Optimal Bayesian experimental design typically involves maximising an objective function
which is the expectation of some appropriately chosen utility function with respect to
the joint distribution of parameters and responses. This objective function is usually not
available in closed form and the design space can be of high dimensionality.

The \code{acebayes} package uses Approximate Coordinate Exchange (ACE; Overstall & Woods, 2015) 
to maximise the expectation of the utility function. 
It does so by approximating the expected utility using Monte Carlo (MC) 
integration. These MC integration evaluations are then used to fit a Gaussian 
process (GP) emulator which can predict the value of the expected utility for
any design. The coordinate and point exchange algorithms are then used to 
maximise the GP emulator prediction of the expected utility. 
}
\details{
\tabular{ll}{
Package: \tab acebayes\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-08-04\cr
License: \tab GPL-2\cr
}
The most important functions are \code{\link{ace}} and \code{\link{aceglm}}.

The function \code{\link{ace}} implements the ACE algorithm. It has two mandatory arguments: \code{utility} (a function specifying the chosen utility function) and \code{start.d} (the initial design). The function will return the final design from the algorithm, along with some information to assess convergence.

The computational time of \code{\link{ace}} is highly dependent on the computational time of evaluating the user-supplied function \code{utility}. Therefore it is recommended that users take advantage of \code{R} packages such as \code{Rcpp} (Eddelbuettel & Francois, 2011), \code{RcppArmadillo} (Eddelbuettel & Sanderson, 2014), or \code{RcppEigen} (Bates & Eddelbuettel, 2013) to implement this function.

The function \code{\link{aceglm}} is a wrapper function for \code{\link{ace}} which uses the ACE algorithm to find pseudo-Bayesian optimal experimental designs for generalised linear models.

For more details on the underpinning methodology, see Overstall & Woods (2015).
}
\author{
Antony M. Overstall & David C. Woods

Maintainer: Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk}
}
\references{
Bates, D. & Eddelbuettel, D. (2013). Fast and Elegant Numerical Linear Algebra Using the 
RcppEigen Package. \emph{Journal of Statistical Software}, \bold{52(5)}, 1-24. 
\url{http://www.jstatsoft.org/v52/i05/}

Eddelbuettel, D. & Francois, R. (2011). Rcpp: Seamless R and C++ Integration. 
\emph{Journal of Statistical Software}, \bold{40(8)}, 1-18. \url{http://www.jstatsoft.org/v40/i08/}

Eddelbuettel, D. & Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance 
C++ linear algebra. \emph{Computational Statistics and Data Analysis}, \bold{71}, 1054-1063.

Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.
}

\keyword{package}

\examples{
## This example uses aceglm to find a Bayesian D-optimal design for a 
## first order logistic regression model with 6 runs 4 factors (i.e. 5 parameters).
## The priors are those used by Overstall & Woods (2015), i.e. a uniform prior 
## distribution is assumed for each parameter. The design space for each coordinate 
## is [-1, 1].

set.seed(1)
## Set seed for reproducibility.

n<-6
## Specify the sample size (number of runs).

start.d<-matrix(2*randomLHS(n=n,k=4)-1,nrow=n,ncol=4,
dimnames=list(as.character(1:n),c("x1","x2","x3","x4")))
## Generate an initial design of appropriate dimension. The initial design is a 
## Latin hypercube sample.

low<-c(-3, 4, 5, -6, -2.5)
upp<-c(3, 10, 11, 0, 3.5)
## Lower and upper limits of the uniform prior distributions.

prior<-function(B){
t(t(6*matrix(runif(n=5*B),ncol=5))+low)}
## Create a function which specifies the prior. This function will return a 
## B by 5 matrix where each row gives a value generated from the prior 
## distribution for the model parameters.

example<-aceglm(formula=~x1+x2+x3+x4, start.d = start.d, family = binomial, 
prior = prior , criterion = "D", B = c(1000,1000), N1 = 1, N2 = 0, upper = 1)
## Call the aceglm function which implements the ACE algorithm requesting 
## only one iteration of Phase I and zero iterations of Phase II. The Monte
## Carlo sample size for the comparison procedure (B[1]) is set to 1000.

example
## Print out a short summary. 

#Generalised Linear Model
#Criterion = Bayesian D-optimality 
#
#Number of runs = 6
#
#Number of factors = 4
#
#Number of Phase I iterations = 1
#
#Number of Phase II iterations = 0
#
#Computer time = 00:00:02

example$phase2.d

#          x1          x2          x3         x4
#1 -0.3571245  0.16069337 -0.61325375  0.9276443
#2 -0.9167309  0.91411512  0.69842151  0.2605092
#3 -0.8843699  0.42863930 -1.00000000 -0.9679402
#4  0.3696224 -0.27126080  0.65284076  0.1850767
#5  0.7172267 -0.34743402 -0.05968457 -0.6588896
#6  0.7469636  0.05854029  1.00000000 -0.1742566
}
