\name{acenlm}
\alias{acenlm}

\title{
Approximate Coordinate Exchange (ACE) Algorithm for Non-Linear Models
}
\encoding{UTF-8}
\description{
This function implements the approximate coordinate exchange algorithm (Overstall & Woods, 2017) for finding optimal Bayesian designs for non-linear models.
}
\usage{
acenlm(formula, start.d, prior, B, criterion = c("D", "A", "E", "SIG", "NSEL"), 
method = c("quadrature", "MC"),  Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, 
progress = FALSE, limits = NULL, tolerence = 1e-08)
}

\arguments{
  \item{formula}{
An object of class \code{"formula"}: a symbolic description of the model. The terms should correspond to the column names of the argument \code{start.d} and the column names of the result of a call to argument \code{prior}.
}
  \item{start.d}{
An \code{n} by \code{k} matrix with column names used by the argument \code{formula}.
}
	\item{prior}{
For \code{method = "MC"}, a function with one argument: \code{B}; a scalar integer. This function should return a \code{B} by \code{p} (or \code{p+1} for argument \code{criterion = "SIG"} and \code{criterion = "NSEL"}) matrix containing a random sample from the prior distribution of the parameters. \code{p} is the number of parameters in \code{formula}. The columns should have names corresponding to the terms in the argument \code{formula}. For \code{criterion = "SIG"} and \code{criterion = "NSEL"}, an extra column (named \code{"sig2"}) should contain a random sample from the prior distribution of the error variance. 

For \code{method = "quadrature"}, a list specifying a normal or uniform prior for the model parameters in \code{formula}. 

For a normal prior distribution, the list should have named entries \code{mu} and \code{sigma2} specifying the prior mean and variance-covariance matrix. The prior mean \code{mu} may be specified as a vector of length \code{p} where the elements are named to match the parameters in the \code{formula} argument. The prior variance-covariance matrix may be specified as either a scalar common variance or a vector of variances (for independent prior distributions) or as a \code{p x p} matrix. The names attribute of \code{mu} must match the names of the terms in the argument \code{formula}.

For a uniform prior distribution, the list should have a named entry \code{limits}, a \code{2 x p} matrix with each column giving the lower and upper limits of the independent continuous uniform distribution for the corresponding parameter. The column names of \code{limits} must match the names of the terms in the argument \code{formula}.

}
 \item{B}{
An optional vector argument of length two. 

For \code{criterion = "MC"}, it specifies the size of the Monte Carlo samples, generated from the joint distribution of model parameters and unobserved responses, to use in the ACE algorithm. The 1st sample size, \code{B[1]}, gives the sample size to use in the comparison procedures, and the 2nd sample size, \code{B[2]}, gives the sample size to use for the evaluations of Monte Carlo integration that are used to fit the Gaussian process emulator. If left unspecfified, the default value is \code{c(20000,1000)}.

For \code{criterion = "quadrature"}, it specifies the tuning parameters (numbers of radial abscissas and random rotations) for the implemented quadrature method; see \bold{Details} for more information. If left unspecified, the default value is \code{c(2, 8)}.

}

	\item{criterion}{
An optional character argument specifying the utility function. There are currently five utility functions implemented consisting of 
\enumerate{
\item
pseudo-Bayesian D-optimality (\code{criterion = "D"}); 
\item
pseudo-Bayesian A-optimality (\code{criterion = "A"}); 
\item
pseudo-Bayesian E-optimality (\code{criterion = "E"}); 
\item
Shannon information gain (\code{criterion = "SIG"}); 
\item
negative squared error loss (\code{criterion = "NSEL"}). 
}
The default value is \code{"D"} denoting pseudo-Bayesian D-optimality. See 'Details' for more information. 
}
	\item{method}{
An optional character argument specifying the method of approximating the expected utility function. Current choices are \code{method = "quadrature"} for a deterministic quadrature approximation and \code{method = "MC"} for a stochastic Monte Carlo approximation. The first of these choices is available for the argument \code{criterion = "A"}, \code{"D"} or \code{"E"}; the second choice is available for all approximations. The default value for A, D, E criteria is \code{"quadrature"}; for SIG and NSL criteria it is MC. See \bold{Details} for more information. 		
}
  \item{Q}{
An integer specifying the number of evaluations of the approximate expected utility (with sample size \code{B[2]}) that are used to fit the Gaussian process emulator. The default value is \code{20}.
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
  \item{limits}{
An argument specifying the grid over which to maximise the Gaussian process emulator for the expected utility function. It should be a function with three arguments: \code{i}, \code{j} and \code{d} which generates a one-dimensional grid for the \code{ij}th coordinate of the design when the current design is \code{d}. The default value is NULL which generates values uniformly on the interval \code{(lower[i,j],upper[i,j])} or \code{(lower,upper)} depending on whether the arguments \code{lower} and \code{upper} are matrices or scalars, respectively.
}
	\item{tolerence}{
For a deterministic approximation, the minimum improvement in approximate expected utility required to proceeed to the next iteration of coordinate (Phase I) or point (Phase II) exchange. The default value is \code{1e-08}. 				
}
}
\details{
This function implements the ACE algorithm to find designs for general classes of nonlinear regression models with identically and independently normally distributed errors. This means the user does not have to write their own utility function.

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
where \eqn{e()}{e()} denotes the function that calculates the eigenvalues of its argument, i.e. the optimal design maximises the prior expectation of the minimum eigenvalue of the Fisher information matrix.}

Two fully Bayesian utility functions are implemented:

\enumerate{
\item
\bold{Shannon information gain} (\code{criterion = "SIG"})
Here the objective function is
\deqn{U^{SIG}(d) = E \left( \pi(\theta|y,d) - \pi(\theta) \right)}{
U^SIG(d) = E ( f(\theta|y,d) - f(\theta) )}
\item
\bold{Negative squared error loss} (\code{criterion = "NSEL"})
Here the objective function is
\deqn{U^{NSEL}(d) = E \left( \left(\theta - E(\theta,y,d)\right)^T \left(\theta - E(\theta,y,d)\right)\right)}{
U^NSEL(d) = E ( (\theta - E(\theta|y,d))^T(\theta - E(\theta|y,d)))}
}

The expected utilities can be approximated using Monte Carlo methods (\code{method = "MC"} for all criteria) or using a deterministic quadrature method (\code{method = "quadrature"}, implemented for the D, A and E criteria). The former approach approximates the expected utility via sampling from the prior. The latter approach uses a radial-spherical integration rule (Monahan and Genz, 1997) and \code{B[1]} specifies the number, \eqn{n_r}{nr}, of radial abscissas and \code{B[2]} specifies the number, \eqn{n_q}{nq}, of random rotations. Larger values of \eqn{n_r}{nr} and \eqn{n_q}{nq} will produce more accurate, but also more computationally expensive, approximations. See Gotwalt et al. (2009) for further details.   


For more details on the ACE algorithm, see Overstall & Woods (2017).
}

\value{
The function will return an object of class \code{"ace"} which is a list with the following components:
\item{utility}{
The utility function resulting from the choice of arguments.
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
A vector containing the evaluations of the approximate expected utility of the current design at each stage of Phase I of the ACE algorithm. This can be used to assess convergence.
}
\item{phase2.trace}{
A vector containing the evaluations of the approximate expected utility of the current design at each stage of Phase II of the ACE algorithm. This can be used to assess convergence.
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
This will be \code{FALSE}.
}
\item{nlm}{
If the object is a result of a direct call to \code{acenlm} then this is \code{TRUE}.
}
\item{criterion}{
If the object is a result of a direct call to \code{acenlm} then this is the argument \code{criterion}.
}
\item{method}{
If the object is a result of a direct call to \code{acenlm} then this is the argument \code{method}.	
}
\item{prior}{
If the object is a result of a direct call to \code{aceglm} then this is the argument \code{prior}.
}
\item{time}{
Computational time (in seconds) to run the ACE algorithm.
}
\item{binary}{
The argument \code{binary}. Will be \code{FALSE} for the utility functions currently implemented.
}
}

\references{
Chaloner, K. & Verdinelli, I. (1995). Bayesian experimental design: a review. \emph{Statistical
Science}, \bold{10}, 273-304.

Gotwalt, C. M., Jones, B. A.  & Steinberg, D. M. (2009). Fast computation of designs robust to parameter uncertainty for nonlinear settings. \emph{Technometrics}, \bold{51}, 88-95.

Meyer, R. & Nachtsheim, C. (1995). The coordinate exchange algorithm for constructing exact
optimal experimental designs. \emph{Technometrics}, \bold{37}, 60-69.

Monahan, J. and Genz, A. (1997). Spherical-radial integration rules for
Bayesian computation,” Journal of the American Statistical Association, 92, 664–674.

Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, in press. (\emph{arXiv:1501.00264}).
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods & Maria Adamou
}

\note{
This is a wrapper function for \code{\link{ace}}.
}

\seealso{
\code{\link{ace}}, \code{\link{aceglm}}.
}
\examples{
## This example uses aceglm to find a Bayesian D-optimal design for a 
## compartmental model with 6 runs 1 factor. The priors are 
## those used by Overstall & Woods (2017). The design space for each 
## coordinate is [0, 24] hours.

set.seed(1)
## Set seed for reproducibility.

n<-6
## Specify the sample size (number of runs).

start.d<-matrix(24*randomLHS(n=n,k=1),nrow=n,ncol=1,
dimnames=list(as.character(1:n),c("t")))
## Generate an initial design of appropriate dimension. The initial design is a 
## Latin hypercube sample.

low<-c(0.01884,0.298,21.8)
upp<-c(0.09884,8.298,21.8)
## Lower and upper limits of the uniform prior distributions. Note that the prior
## for the third element is a point mass.

prior<-function(B){
out<-cbind(runif(n=B,min=low[1],max=upp[1]),runif(n=B,min=low[2],max=upp[2]),
runif(n=B,min=low[3],max=upp[3]))
colnames(out)<-c("a","b","c")
out}

## Create a function which specifies the prior. This function will return a 
## B by 3 matrix where each row gives a value generated from the prior 
## distribution for the model parameters.

example1<-acenlm(formula=~c*(exp(-a*t)-exp(-b*t)), start.d = start.d, prior = prior, 
N1 = 1, N2 = 0, B = c(1000,1000), lower = 0, upper = 24, method = "MC")
## Call the acenlm function which implements the ACE algorithm requesting 
## only one iteration of Phase I and zero iterations of Phase II. The Monte
## Carlo sample size for the comparison procedure (B[1]) is set to 1000.

example1
## Print out a short summary.

#Non Linear Model 
#Criterion = Bayesian D-optimality 

#Number of runs = 6

#Number of factors = 1

#Number of Phase I iterations = 1

#Number of Phase II iterations = 0

#Computer time = 00:00:00

example1$phase2.d
## Look at the final design.

#           t
#1 19.7787011
#2  2.6431912
#3  0.2356938
#4  8.2471451
#5  1.4742319
#6 12.7062270

prior2 <- list(limits = cbind(rbind(low, upp)))
colnames(prior2$limits) <- c("a", "b", "c")
example2 <- acenlm(formula = ~c * (exp(-a * t) - exp(-b *t)), start.d = start.d, prior = prior2, 
lower = 0, upper = 24, N1 = 1, N2 = 0)
## Call the acenlm function with the default method of "quadrature"

example2$phase2.d
## Final design

#           t
#1  0.3090187
#2  2.0710214
#3  6.5164562
#4  0.5964350
#5 23.3100644
#6 20.1339899

utility <- utilitynlm(formula = ~c * (exp(-a * t) - exp(-b *t)), prior = prior, 
                            desvars = "t", method = "MC")$utility
## create a utility function to compare designs

mean(utility(example1$phase1.d, B = 20000))
#[1] 12.13958
mean(utility(example2$phase1.d, B = 20000))
#[1] 12.10159
## Compare the two designs using the Monte Carlo approximation
}
