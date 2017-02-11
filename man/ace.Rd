\name{ace}
\alias{ace}
\alias{acephase1}
\alias{acephase2}

\title{
Approximate Coordinate Exchange (ACE) Algorithm
}
\encoding{UTF-8}
\description{
This function implements the approximate coordinate exchange (ACE) algorithm (Overstall & Woods, 2017) for finding optimal Bayesian experimental designs by maximising an approximation to an intractable expected utility function.
}
\usage{
ace(utility, start.d, B, Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, limits = NULL, 
progress = FALSE, binary = FALSE, deterministic = FALSE, tolerence = 1e-08)
		
acephase1(utility, start.d, B, Q = 20, N1 = 20, lower, upper, limits = NULL, 
progress = FALSE, binary = FALSE, deterministic = FALSE, tolerence = 1e-08)

acephase2(utility, start.d, B, N2 = 100, progress = FALSE, binary = FALSE, 
deterministic = FALSE, tolerence = 1e-08)
}

\arguments{
  \item{utility}{
A function with two arguments: \code{d} and \code{B}. For a Monte Carlo approximation, it should return a \emph{\bold{vector}} of length \code{B} where each element gives the value of the utility function for design \code{d}, for a value generated from the joint distribution of model parameters and unobserved responses. For a deterministic (e.g. quadrature) approximation, it should return a \emph{\bold{scalar}} giving the approximate value of the expected utility for design \code{d}. In this latter case, the argument \code{B} can be a list containing tuning parameters for the deterministic approximation (see the examples). If \code{B} is not required, the utility function must still accept the argument, e.g. using the \code{...} notation.
}
  \item{start.d}{
An \code{n} by \code{k} matrix specifying the initial design for the ACE algorithm.
}
  \item{B}{
For a Monte Carlo approximation, a vector of length two specifying the size of the Monte Carlo samples, generated from the joint distribution of model parameters and unobserved responses, to use in the ACE algorithm. The 1st sample size, \code{B[1]}, gives the sample size to use in the comparison procedures, and the 2nd sample size, \code{B[2]}, gives the sample size to use for the evaluations of Monte Carlo integration that are used to fit the Gaussian process emulator. If missing when \code{deterministic = FALSE}, the default value is \code{c(20000,1000)}. For a deterministic approximation, \code{B} may be a list of length two containing any necessary tuning parameters for the utility calculations for the comparison and emulation steps (see examples).
}
  \item{Q}{
An integer specifying the number of evaluations of the approximate expected utility (e.g. Monte Carlo integration with sample size \code{B[2]}) that are used to fit the Gaussian process emulator. The default value is \code{20}.
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
An argument specifying the bounds on the design space. This argument can either be a scalar or a matrix of the same dimension as the argument \code{start.d} which specifies the lower limits of all coordinates of the design space. The default value is -1.
}
  \item{upper}{
An argument specifying the bounds on the design space. This argument can either be a scalar or a matrix of the same dimension as the argument \code{start.d} which specifies the upper limits of all coordinates of the design space. The default value is 1.
}
  \item{limits}{
An argument specifying the grid over which to maximise the Gaussian process emulator for the expected utility function. It should be a function with three arguments: \code{i}, \code{j} and \code{d} which generates a one-dimensional grid for the \code{ij}th coordinate of the design when the current design is \code{d}. The default value is NULL which generates values uniformly on the interval \code{(lower[i,j],upper[i,j])} or \code{(lower,upper)} depending on whether the arguments \code{lower} and \code{upper} are matrices or scalars, respectively.
}
  \item{progress}{
A logical argument indicating whether the iteration number and other information detailing the progress of the algorithm should be printed. The default value is \code{FALSE}.
}
  \item{binary}{
A logical argument indicating whether the utility function has binary or continuous output. In some cases, the utility function is an indicator function of some event giving binary output. The expected utility function will then be the expected posterior probability of the event. Utility functions such as Shannon information gain and negative squared error loss give continuous output. The type of output guides the choice of comparison procedure used in the ACE algorithm. The default value is \code{FALSE}, indicating the utility function has continuous output.
}
	\item{deterministic}{
A logical argument indicating if a Monte Carlo (\code{FALSE}, default) or deterministic (\code{TRUE}) approximation to the expected utility is being used. 		
}
	\item{tolerence}{
For a deterministic approximation, the minimum improvement in approximate expected utility required to proceeed to the next iteration of coordinate (Phase I) or point (Phase II) exchange. The default value is \code{1e-08}. 		
}
}
\details{
Optimal Bayesian experimental design (Chaloner & Verdinelli, 1995) involves maximising 
an objective function given by the expectation of some appropriately chosen utility 
function with respect to the joint distribution of parameters and responses, over a 
design space of all possible designs. This objective function is usually not available 
in closed form and the design space can be of high dimensionality.

Overstall & Woods (2017) proposed the approximate coordinate exchange (ACE) algorithm 
to approximately maximise the expectation of the utility function. The ACE consists of
two phases.

Phase I uses a continuous version of the coordinate exchange algorithm (Meyer & 
Nachtsheim, 1995) to maximise an approximation to the expected loss. Very briefly, 
the approximate expected loss is sequentially minimised over each one-dimensional element
of the design space. The approximate expected loss is given by the predictive mean of a 
Gaussian process (GP) regression model (also known as an emulator or surrogate) fitted 
to a 'small' number (argument \code{Q}) of evaluations of either a Monte Carlo (MC) or deterministic (e.g. quadrature) approximation to the expected utility (the MC sample size or arguments for the deterministic approximation are given by \code{B}). A GP 
emulator is a statistical model and, similar to all statistical models, can be an 
inadequate representation of the underlying process (i.e. the expected utility). 
Instead of automatically accepting the new design given by the value that maximises 
the GP emulator, for MC approximations a Bayesian hypothesis test, independent of the GP emulator, is 
performed to assess whether the expected utility of the new design is larger than the 
current design. For deterministic approximations, the approximate expected utility is calculated for the new design, and compared to that for the current design.

Phase I tends to produce clusters of design points. This is where, for example, two
design points are separated by small Euclidean distance. Phase II allows these points 
to be consolidated into a single repeated design point by using a point exchange 
algorithm (e.g Gotwalt et al., 2009) with a candidate set given by the final design 
from Phase I. In the same way as Phase I, comparisons of the expected loss between 
two designs is made on the basis of either a Bayesian hypothesis test or a direct comparison of deterministic approximations.

The original Bayesian hypothesis test proposed by Overstall & Woods (2017) is 
appropriate for utility functions with continuous domain. Overstall et al. (2017) 
extended the idea to utility functions with binary domain, e.g. the utility 
function is an indicator function for some event.
}

\value{
The function will return an object of class \code{"ace"} which is a list with the following components:

\item{utility}{
The argument \code{utility}.
}
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
A vector containing the approximated expected utility of the current design at each stage of Phase I of the ACE algorithm. This can be used to assess convergence for MC approximations.
}
\item{phase2.trace}{
A vector containing the approximated expected utility of the current design at each stage of Phase II of the ACE algorithm. This can be used to assess convergence for MC approximations.
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
\item{nlm}{
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
\item{deterministic}{
The argument \code{deterministic}.	
}
\item{tolerence}{
The argument \code{tolerence}.	
}
}
\references{
Chaloner, K. & Verdinelli, I. (1995). Bayesian Experimental Design: A Review. \emph{Statistical
Science}, \bold{10}, 273-304.

Gotwalt, C., Jones, B. & Steinberg, D. (2009). Fast Computation of Designs Robust to Parameter
Uncertainty for Nonlinear Settings. \emph{Technometrics}, \bold{51}, 88-95.

Meyer, R. & Nachtsheim, C. (1995). The Coordinate Exchange Algorithm for Constructing Exact
Optimal Experimental Designs. \emph{Technometrics}, \bold{37}, 60-69.

Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, in press. (\emph{arXiv:1501.00264}).

Overstall, A.M., McGree, J.M. & Drovandi, C.C. (2017). An approach for finding fully Bayesian optimal designs using normal-based approximations to loss functions. \emph{Statistics and Computing}, in press, \emph{arXiv:1608.05815}.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods & Maria Adamou
}

\note{
For more detail on the \code{R} implementation of the utility function used in the example below, see \code{\link{utilcomp18bad}}.
}

\examples{
set.seed(1)
## Set seed for reproducibility.

## This example involves finding a Bayesian D-optimal design for a 
## compartmental model with n = 18 runs. There are three parameters. 
## Two parameters have uniform priors and the third has a prior 
## point mass. For more details see Overstall & Woods (2017).

start.d<-optimumLHS(n=18,k=1)
## Create an initial design.

## Using a MC approximation
example1<-ace(utility = utilcomp18bad, start.d = start.d, N1 = 1, N2 = 2, B = c(100,20))
## Implement the ACE algorithm with 1 Phase I iterations and 2 Phase II
## iterations. The Monte Carlo sample sizes are 100 (for comparison) and 20 for
## fitting the GP emulator.

example1
## Produce a short summary.

#User-defined model & utility 
#
#Number of runs = 18
#
#Number of factors = 1
#
#Number of Phase I iterations = 1
#
#Number of Phase II iterations = 2
#
#Computer time = 00:00:00

mean(utilcomp18bad(d=example1$phase2.d,B=100))
## Calculate an approximation to the expected utility for the final design.
## Should get:

#[1] 10.34614

\dontrun{
plot(example1)
## Produces a trace plot of the current value of the expected utility. This
## can be used to assess convergence.
}

## Using a deterministic approximation with the log determinant of the Fisher information matrix 
## averaged across a fixed Latin hypercube design with 25 points
Q <- randomLHS(25, 2)
Q[, 1] <- 0.01884 + 0.298 * Q[, 1]
Q[, 2] <- 0.09884 + 0.298 * Q[, 2]
utilcomp18badF <- function(d, B) {
	mean(as.vector(.Call("utilcomp18badcpp", d, B, PACKAGE = "acebayes")))	
}
example2 <- ace(utility = utilcomp18badF, start.d = start.d, N1 = 1, N2 = 2, B = list(Q, Q), 
deterministic = TRUE)
## Here, B = list(Q, Q) means the same deterministic approximation is used for both the 
## comparison and emulation steps

example2
#User-defined model & utility 
#
#Number of runs = 18
#
#Number of factors = 1
#
#Number of Phase I iterations = 1
#
#Number of Phase II iterations = 2
#
#Computer time = 00:00:00

mean(utilcomp18bad(d=example2$phase2.d,B=100))
#[1] 0.7127533
}