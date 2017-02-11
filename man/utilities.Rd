\name{utilities}
\alias{utilitynlm}
\alias{utilityglm}

\title{
Approximate expected utility function for generalised linear models and non-linear regression models
}
\encoding{UTF-8}
\description{
Generates an approximate utility function for generalised linear models and non-linear regression models.
}
\usage{
utilityglm(formula, family, prior, 
criterion = c("D", "A", "E", "SIG", "NSEL", "SIG-Norm", "NSEL-Norm"),
method = c("quadrature", "MC"), nrq)

utilitynlm(formula, prior, desvars, criterion =  c("D", "A", "E", "SIG", "NSEL"),
method = c("quadrature", "MC"), nrq)
}

\arguments{
  \item{formula}{
An object of class \code{"formula"}: a symbolic description of the model. The terms should correspond to the column names of the argument \code{start.d} and the column names of the result of a call to argument \code{prior}.
}
	\item{family}{
For \code{utilityglm}, a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
}
	\item{prior}{
For \code{method = "MC"}, a function with one argument: \code{B}; a scalar integer. This function should return a \code{B} by \code{p} (or \code{p+1} for when the argument \code{criterion = "SIG"} or \code{criterion = "NSEL"}) matrix containing a random sample from the prior distribution of the parameters. \code{p} is the number of parameters in \code{formula}. The columns should have names corresponding to the terms in the argument \code{formula}. For \code{criterion = "SIG"} and \code{criterion = "NSEL"}, an extra column (named \code{"sig2"}) should contain a random sample from the prior distribution of the error variance. 

For \code{method = "quadrature"}, a list specifying a normal or uniform prior for the model parameters in \code{formula}. 

For a normal prior distribution, the list should have named entries \code{mu} and \code{sigma2} specifying the prior mean and variance-covariance matrix. The prior mean may be specified as a vector of length \code{p}. The prior variance-covariance matrix may be specified as either a scalar common variance or a vector of variances (for independent prior distributions) or as a \code{p x p} matrix. The names attribute of \code{mu} must match the names of the terms in the argument \code{formula}.

For a uniform prior distribution, the list should have a named entry \code{limits}, a \code{2 x p} matrix with each column giving the lower and upper limits of the independent continuous uniform distribution for the corresponding parameter. The column names of \code{limits} must match the names of the terms in the argument \code{formula}.

}

\item{desvars}{
For \code{utilitynlm}, a character vector listing the design variables that appear in the argument \code{formula}.
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
An optional character argument specifying the method of approximating the expected utility function. Current choices are \code{method = "quadrature"} for a deterministic quadrature approximation and \code{method = "MC"} for a stochastic Monte Carlo approximation. The first of these choices is available for argument \code{criterion = "A"}, \code{"D"} or \code{"E"}; the second choice is available for all approximations. The default value for A, D, E criteria is \code{"quadrature"}; for SIG and NSL criteria it is MC. See \bold{Details} for more information. 		
}
\item{nrq}{
For \code{method = "quadrature"}, a vector of length two specifying the number of radial abscissas (\code{nrq[1]}) and random rotations (\code{nrq[2]}) required for the implemented quadrature scheme; see \bold{Details} for more information. If left unspecified, the default value is \code{c(2, 8)}.	
}
}
\details{
A function that generates approximate expected utility functions for Bayesian design (Chaloner & Verdinelli, 1995).
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
Note that these two fully Bayesian utility functions are only implemented for logistic regression, i.e. \code{family = binomial(link="logit")} and Poisson regression, i.e. \code{family = poisson(link="log")}.

The expected utilities can be approximated using Monte Carlo methods (\code{method = "MC"} for all criteria, as used by Overstall & Woods (2017) or using a deterministic quadrature method (\code{method = "quadrature"}, implemented for the D, A and E criteria). The former approach approximates the expected utility via sampling from the prior. The latter approach uses a radial-spherical integration rule (Monahan and Genz, 1997) and \code{B[1]} specifies the number, nr, of radial abscissas and \code{B[2]} specifies the number, nq, of random rotations. Larger values of nr and nq will produce more accurate, but also more computationally expensive, approximations. See Gotwalt et al. (2009) for further details.   


}

\value{
The function will return a list with the following components:
\item{utility}{
The utility function resulting from the choice of arguments.
}
}

\references{
Chaloner, K. & Verdinelli, I. (1995). Bayesian experimental design: a review. \emph{Statistical
Science}, \bold{10}, 273-304.

Gotwalt, C. M., Jones, B. A.  & Steinberg, D. M. (2009). Fast computation of designs robust to parameter uncertainty for nonlinear settings. \emph{Technometrics}, \bold{51}, 88-95.

Monahan, J. and Genz, A. (1997). Spherical-radial integration rules for
Bayesian computation,” Journal of the American Statistical Association, 92, 664-674.

Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, in press. (\emph{arXiv:1501.00264}).
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods & Maria Adamou
}

\seealso{
\code{\link{aceglm}},\code{\link{acenlm}}.
}
\examples{
## 1. This example uses utilityglm to generate the pseudo-Bayesian D-optimality
## approximate expected utility function using a Monte Carlo approximation.

low<-c(-3, 4, 5, -6, -2.5)
upp<-c(3, 10, 11, 0, 3.5)
## Lower and upper limits of the uniform prior distributions.

prior<-function(B){
t(t(6*matrix(runif(n=5*B),ncol=5))+low)}
## Create a function which specifies the prior. This function will return a 
## B by 5 matrix where each row gives a value generated from the prior 
## distribution for the model parameters.

ex <- utilityglm(formula = ~x1+x2+x3+x4, family = binomial, prior = prior, method = "MC")

set.seed(1)
## Set seed for reproducibility.

n<-6
## Specify the sample size (number of runs).

start.d<-matrix(2*randomLHS(n=n,k=4)-1,nrow=n,ncol=4,
dimnames=list(as.character(1:n),c("x1","x2","x3","x4")))
## Generate an initial design of appropriate dimension. The initial design is a 
## Latin hypercube sample.

ex$utility(d = start.d, B = 10)
## Evaluate resulting approximate utility. Should get:

## [1] -14.01218 -18.53715 -20.78817 -24.02731 -14.89368 -14.25400 -17.38152
## [8] -17.70984 -14.15254 -21.40505

## 2. This example uses utilitynlm to generate the psuedo-Bayesian A-optimality expected utility
## function using a quadrature approximation

low<-c(0.01884,0.298,21.8)
upp<-c(0.09884,8.298,21.8)
## Lower and upper limits of the uniform prior distributions. Note that the prior
## for the third element is a point mass.

prior2 <- list(limits = cbind(rbind(low, upp)))
colnames(prior2$limits) <- c("a", "b", "c")
## Specify a uniform prior with ranges given by low and upp

ex2 <- utilitynlm(formula = ~c * (exp(-a * t) - exp(-b *t)), prior = prior2, 
                            desvars = "t")
                            
n <- 6
start.d <- matrix(24 * randomLHS(n = n, k = 1), nrow = n)
colnames(start.d) <- "t"
ex2$utility(d = start.d) 
## -3.651616    

}
