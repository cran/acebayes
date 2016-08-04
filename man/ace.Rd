\name{ace}
\alias{ace}

\title{
Approximate Coordinate Exchange (ACE) Algorithm
}
\encoding{UTF-8}
\description{
This function implements the approximate coordinate exchange algorithm (Overstall & Woods, 2015) for finding optimal Bayesian experimental designs by approximately maximising an intractable expected utility function.
}
\usage{
ace(utility, start.d, B = c(20000,1000), Q = 20, 
		N1 = 20, N2 = 100, lower = -1, upper = 1, 
		limits = NULL, progress = FALSE, binary = FALSE)
}

\arguments{
  \item{utility}{
A function with two arguments: \code{d} and \code{B}. It should return a vector of length \code{B} where each element gives the value of the utility function for design \code{d}, for a value generated from the joint distribution of model parameters and unobserved responses.
}
  \item{start.d}{
An \code{n} by \code{k} matrix specifying the initial design for the ACE algorithm.
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
  \item{limits}{
An argument specifying the grid over which to maximise the Gaussian process emulator for the expected utility function. It should be a function with three arguments: \code{i}, \code{j} and \code{d} which generates a one-dimensional grid for the \code{ij}th coordinate of the design when the current design is \code{d}. The default value is NULL which generates values uniformly on the interval \code{(lower[i,j],upper[i,j])}.
}
  \item{progress}{
A logical argument indicating whether the iteration number should be printed. The default value is \code{FALSE}.
}
  \item{binary}{
A logical argument indicating whether the utility function outputs binary or real numbers. In some cases, the utility function is an indicator function of some event. The expected utility function will then be the expected posterior probability of the event. Utility functions such as Shannon information gain and negative squared error loss output real numbers. The default value is \code{FALSE} - indicating the utility function outputs real numbers.
}
}
\details{
Optimal Bayesian experimental design typically involves maximising an objective function
which is the expectation of some appropriately chosen utility function with respect to
the joint distribution of parameters and responses. This objective function is usually not
available in closed form and the design space can be of high dimensionality.

Overstall & Woods (2015) proposed the approximate coordinate exchange algorithm (ACE) to 
approximately maximise the objective function. The algorithm has two phases. Phase I involves 
using the coordinate exchange algorithm (Meyer & Nachtsheim, 1995) to maximise the objective 
function but where a Gaussian process (GP) emulator has been used to approximate the objective 
function. The GP emulator is based on a small number of Monte Carlo integration evaluations of 
the objective function. Phase II involves using the point exchange algorithm (see, for example, 
Morris 2011, pages 310-311) to consolidate clusters of design points.

See Overstall & Woods (2015) for more details.
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
A vector containing the Monte Carlo integration evaluations of the expected utility of the current design at each stage of Phase II of the ACE algorithm. This can be used to assess convergence.
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
If the object is a result of a direct call to \code{ace} then this is \code{FALSE}.
}
\item{criterion}{
If the object is a result of a direct call to \code{ace} then this is \code{"NA"}.
}
\item{prior}{
If the object is a result of a direct call to \code{ace} then this is \code{"NA"}.
}
\item{time}{
Computational time (in seconds) to run the ACE algorithm.
}
\item{binary}{
The argument \code{binary}.
}
}
\references{
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
For more detail on the \code{R} implementation of the utility function used in the example below, see \code{\link{utilcomp18bad}}.
}

\seealso{
\code{\link{acephase1}},
\code{\link{acephase2}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility.

## This example involves finding a Bayesian D-optimal design for a 
## compartmental model with n = 18 runs. There are three parameters. 
## Two parameters have uniform priors and the third has a prior 
## point mass. For more details see Overstall & Woods (2015).

start.d<-optimumLHS(n=18,k=1)
## Create an initial design.

example<-ace(utility = utilcomp18bad, start.d = start.d, N1 = 1, N2 = 2, B = c(100,20))
## Implement the ACE algorithm with 1 Phase I iterations and 2 Phase II
## iterations. The Monte Carlo sample sizes are 100 (for comparison) and 20 for
## fitting the GP emulator.

example
## Produce a short summary.

#User-defined utility 
#
#Number of runs = 18
#
#Number of factors = 1
#
#Number of Phase I iterations = 1
#
#Number of Phase II iterations = 2
#
#Computer time = 00:00:01

mean(utilcomp18bad(d=example$phase2.d,B=100))
## Calculate an approximation to the expected utility for the final design.
## Should get:

#[1] 10.34614

\dontrun{
plot(example)
## Produces a trace plot of the current value of the expected utility. This
## can be used to assess convergence.
}
}