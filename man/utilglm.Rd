\name{utilglm}
\alias{utilglm}

\title{
Pseudo-Bayesian Utility Functions for Generalised Linear Models
}
\encoding{UTF-8}
\description{
This function calculates the pseudo-Bayesian utility functions for generalised linear models. 
In particular, the utility functions that result in Bayesian D-optimality, 
Bayesian A-optimality, and Bayesian E-optimality.
}
\usage{
utilglm(x, beta, family, criterion = "D")
}

\arguments{
  \item{x}{
An n by k model matrix.
}
  \item{beta}{
A B by k matrix where each row contains a value generated from the prior distribution of the
model parameters.
}
	\item{family}{
A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
}
	\item{criterion}{
An optional character specifiying the Bayesian optimality utility function. See 'Details' for more information.
}
}	
\details{
Chaloner & Verdinelli (1995) describe the relationship between approximations to the expectation of commonly-used fully Bayesian utility functions and classical optimality criteria for non-linear models (which produce designs sometimes known as pseudo-Bayesian optimal designs).

Typically, the objective function is the prior expectation of some scalar summary of the Fisher information matrix. This objective function is usually intractable. This function implements the ACE algorithm to approximately maximise these objective functions for GLMs thus finding pseudo-Bayesian optimal designs.

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
The functions will return a vector of length B, where each element gives the value of the utility function for the 
model parameters in the corresponding row of \code{beta}.
}

\references{
Chaloner, K. & Verdinelli, I. (1995). Bayesian Experimental Design: A Review. \emph{Statistical
Science}, \bold{10}, 273-304.

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
\code{\link{aceglm}}
}
\examples{
set.seed(1)
## Set seed for reproducibility

x <- cbind( 1 , runif(10) )
## Create a model matrix

beta <- cbind( runif(10) , runif(10) )
## Generate a sample of size 10 from a prior distribution

utilglm(x = x, beta = beta, family = binomial(link="cloglog"), criterion = "D")
## For logistic regression, evaluate utilty function for corresponding to 
## Bayesian D-optimality for each of the 10 values in the prior sample. 

# [1] 1.0877268 1.2850148 0.7216932 1.3221948 1.0203772 1.2231432 1.2356192
# [8] 0.3998329 0.9629448 0.9380129 
}
