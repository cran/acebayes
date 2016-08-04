\name{aceglm}
\alias{aceglm}

\title{
Approximate Coordinate Exchange (ACE) Algorithm for Generalised Linear Models
}
\encoding{UTF-8}
\description{
This function implements the approximate coordinate exchange algorithm (Overstall & Woods, 2015) for finding optimal pseudo-Bayesian designs for generalised linear models (GLMs).
}
\usage{
aceglm(formula, start.d, family, prior, criterion = "D", 
B = c(20000,1000),  Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, 
progress = FALSE)
}

\arguments{
  \item{formula}{
An object of class \code{"formula"}: a symbolic description of the model. The terms should correspond to the column names of the argument \code{start.d}.
}
  \item{start.d}{
An \code{n} by \code{k} matrix with column names used by the argument \code{formula}.
}
	\item{family}{
A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
}
	\item{prior}{
A function with one argument: \code{B}; a scalar integer. This function should return a \code{B} by \code{p} matrix containg a random sample from the prior distribution of the parameters. 
}
	\item{criterion}{
An optional character argument specifying the pseudo-Bayesian utility function. See 'Details' for more information. The default value is \code{"D"}.
}
  \item{B}{
A vector of length two specifying the size of the samples, generated from the joint distribution of model parameters and unobserved responses, to use in the ACE algorithm. The 1st sample size, \code{B[1]}, gives the sample size to use in the comparison procedures, and the 2nd sample size, \code{B[2]}, gives the sample size to use for the evaluations of Monte Carlo integration that are used to fit the Gaussian process emulator. The default value is \code{c(20000,1000)}.
}
  \item{Q}{
An integer specifying the number of evaluations of Monte Carlo integration (with sample size \code{B[2]}) that are used to fit the Gaussian process emulator. The default value is \code{20}.
}
  \item{N1}{
An integer specifying the number of iterations of Phase I of the ACE algorithm (the coordinate exchange phase). 
The default value is \code{20}.
}
  \item{N2}{
An integer specifying the number of iterations of Phase II of the ACE algorithm (the point exchange phase).
The default value is \code{100}. 
}
  \item{lower}{
An argument specifying the design space. This argument can either be a scalar or a matrix of the same dimension as the argument \code{start.d} which specifies the lower limits of all coordinates of the design space. The default value is -1.
}
  \item{upper}{
An argument specifying the design space. This argument can either be a scalar or a matrix of the same dimension as the argument \code{start.d} which specifies the upper limits of all coordinates of the design space. The default value is 1.
}
  \item{progress}{
 A logical argument indicating whether the iteration number should be printed. The default value is \code{FALSE}.
}
}
\details{
Chaloner & Verdinelli (1995) describe the relationship between approximations to the expectation of commonly-used fully Bayesian utility functions and classical optimality criteria for non-linear models (which produce designs sometimes known as pseudo-Bayesian optimal designs).

Typically, the objective function is the prior expectation of some scalar summary of the Fisher information matrix. This objective function is usually intractable. This function implements the ACE algrorithm to approximately maximise these objective functions for GLMs thus finding pseudo-Bayesian optimal designs.

Three different pseudo-Bayesian utility functions (or criteria) are implemented:

\enumerate{
\item
\bold{Bayesian D-optimal} (\code{criterion = "D"})

Here the objective function is
\deqn{U^D(d) = E \left( \log \vert \mathcal{I} (\theta;d) \vert \right)}{
U^D(d) = E ( log | I(\theta;d)|)}
i.e. the optimal design maximises the prior expectation of the log determinant of the Fisher information matrix.
\item
\bold{Bayesian A-optimal} (\code{criterion = "A"})

Here the objective function is
\deqn{U^A(d) = - E \left( \mathrm{tr} \left\{ \mathcal{I} (\theta;d)^{-1} \right\} \right)}{
U^A(d) = - E ( tr (I(\theta;d)^(-1)))}
i.e. the optimal design minimises the prior expectation of the trace of the inverse Fisher information matrix.
\item
\bold{Bayesian E-optimal} (\code{criterion = "E"})

Here the objective function is
\deqn{U^E(d) = E \left( \mathrm{min} \mbox{ } e\left(\mathcal{I} (\theta;d) \right) \right)}{
U^E(d) = E ( min(e(I(\theta;d))))}
where \eqn{e()}{e()} denotes the function that calculates the eigenvalues of its argument, i.e. the optimal design maximises the prior expectation of the minimum eigenvalue of the Fisher information matrix.
}

For more details on the ACE algorithm and pseudo-Bayesian designs, see Overstall & Woods (2015).
}

\value{
The function will return an object of class \code{"ace"} which is a list with the following components:

\item{start.d}{
The argument \code{start.d}.
}
\item{phase1.d}{
The design found from Phase I of the ACE algorithm.
}
\item{phase2.d}{
The design found from Phase II of the ACE algorithm.
}
\item{phase1.trace}{
A vector containing the Monte Carlo integration evaluations of the expected utility of the current design at each stage of Phase I of the ACE algorithm. This can be used to assess convergence.
}
\item{phase2.trace}{
A vector containing the Monte Carlo inegration evaluations of the expected utility of the current design at each stage of Phase II of the ACE algorithm. This can be used to assess convergence.
}
\item{Q}{
The argument \code{Q}.
}
\item{N1}{
The argument \code{N1}.
}
\item{N2}{
The argument \code{N2}.
}
\item{glm}{
If the object is a result of a direct call to \code{aceglm} then this is \code{TRUE}.
}
\item{criterion}{
If the object is a result of a direct call to \code{aceglm} then this is the argument \code{criterion}.
}
\item{prior}{
If the object is a result of a direct call to \code{aceglm} then this is the argument \code{prior}.
}
\item{time}{
Computational time (in seconds) to run the ACE algorithm.
}
\item{binary}{
The argument \code{binary} to the underlying call to \code{ace}. It indicates whether the utility function outputs binary or real numbers. In some cases, the utility function is an indicator function of some event. The expected utilty function will then be the expected posterior probability of the event. Utility functions such as Shannon information gain and negative squared error loss output real numbers. Since the pseudo-Bayesian utility functions all output real numbers it will be \code{FALSE}.
}
}

\references{
Chaloner, K. & Verdinelli, I. (1995). Bayesian Experimental Design: A Review. \emph{Statistical
Science}, \bold{10}, 273-304.

Meyer, R. & Nachtsheim, C. (1995). The Coordinate Exchange Algorithm for Constructing Exact
Optimal Experimental Designs. \emph{Technometrics}, \bold{37}, 60-69.

Morris, M. (2011). \emph{Design of Experiments: An Introduction Based on Linear Models}. Chapman
and Hall.

Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.
}

\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk} & David C. Woods.
}

\note{
This is a wrapper function for \code{\link{ace}}.
}

\seealso{
\code{\link{ace}}.
}
\examples{
## This example uses aceglm to find a Bayesian D-optimal design for a 
## first order logistic regression model with 6 runs 4 factors. The priors are 
## those used by Overstall & Woods (2015), i.e. each of the five elements of the
## parameters is uniform. The design space for each coordinate is [-1, 1].

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
prior = prior , N1 = 1, N2 = 0, B = c(1000,1000), upper = 1)
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
## Look at the final design.

#          x1          x2          x3         x4
#1 -0.3571245  0.16069337 -0.61325375  0.9276443
#2 -0.9167309  0.91411512  0.69842151  0.2605092
#3 -0.8843699  0.42863930 -1.00000000 -0.9679402
#4  0.3696224 -0.27126080  0.65284076  0.1850767
#5  0.7172267 -0.34743402 -0.05968457 -0.6588896
#6  0.7469636  0.05854029  1.00000000 -0.1742566
}
