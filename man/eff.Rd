\name{eff}
\alias{effnlm}
\alias{effglm}

\title{
Approximate pseudo-Bayesian relative efficiency for generalised linear models and non-linear regression models
}
\encoding{UTF-8}
\description{
Calculates an approximate pseudo-Bayesian relative efficiency of two designs assuming a generalised linear model or non-linear regression model for the response.
}
\usage{
effglm(d1, d2, B, formula=NULL, family=NULL, prior=NULL, 
criterion = c("D", "A", "E"), method = c("quadrature", "MC"))

effnlm(d1, d2, B, formula=NULL, prior=NULL, 
criterion =  c("D", "A", "E"), method = c("quadrature", "MC"))
}

\arguments{
\item{d1, d2}{
Either objects of class \code{"ace"} or \code{n} by \code{k} matrices. The efficiency of \code{d1} will be calculated relative to \code{d2}.
}
\item{B}{
An optional argument for controlling the approximation to the expected utility. 

For \code{criterion = "MC"}, it specifies the size of the Monte Carlo samples, generated from the joint distribution of unknown quantities. If left unspecified, the default value is \code{20000}.

For \code{criterion = "quadrature"}, it is a vector of length \code{2}, specifying the tuning parameters (numbers of radial abscissas and random rotations) for the implemented quadrature method. If left unspecified, the default value is \code{c(2, 8)}.
}
  \item{formula}{
An optional argument providing a symbolic description of the model.

For \code{effglm}, an object of class \code{"formula"}: a symbolic description of the model. The terms should correspond to the column names of the arguments \code{d1} and \code{d2}.

For \code{effnlm}, an object of class \code{"formula"}: a symbolic description of the model. The terms should correspond to the column names of the arguments \code{d1} and \code{d2}, and the argument \code{prior}
}
	\item{family}{
For \code{effglm}, an optional description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
}
	\item{prior}{
An optional argument specifying the prior distribution.  	
	
For \code{method = "MC"}, a function with one argument: \code{B}; a scalar integer. This function should return a \code{B} by \code{p} matrix, with \code{p} the number of model parameters, containing a random sample from the prior distribution of the parameters. The value of \code{p} should correspond to the number of terms specified by the \code{formula} argument. For \code{effnlm}, the column names must match the names of parameters in the \code{formula} argument. 

For \code{method = "quadrature"}, a list specifying a normal or uniform prior for the model parameters. For a normal prior distribution, the list should have named entries \code{mu} and \code{sigma2} specifying the prior mean and variance-covariance matrix. The prior mean may be specified as a scalar, which will then be replicated to form an vector of common prior means, or a vector of length \code{p}. The prior variance-covariance matrix may be specified as either a scalar common variance or a vector of length \code{p} of variances (for independent prior distributions) or as a \code{p} by \code{p} matrix. For \code{effnlm}, the names attribute of \code{mu} must match the names of the parameters in the \code{formula} argument. For a uniform prior distribution, the list should have a named entry \code{support}, a \code{2} by \code{p} matrix with each column giving the lower and upper limits of the support of the independent continuous uniform distribution for the corresponding parameter. For \code{effnlm}, the column names of \code{support} must match the names of parameters in the \code{formula} argument.
}

	\item{criterion}{
An optional character argument specifying the utility function. There are currently three utility functions implemented consisting of 
\enumerate{
\item
\code{pseudo-Bayesian D-optimality} (\code{criterion = "D"}); 
\item
\code{pseudo-Bayesian A-optimality} (\code{criterion = "A"}); 
\item
\code{pseudo-Bayesian E-optimality} (\code{criterion = "E"}).
}
The default value is \code{"D"} denoting pseudo-Bayesian D-optimality. See 'Details' for more information. 
}
\item{method}{
An optional character argument specifying the method of approximating the expected utility function. Current choices are \code{method = "quadrature"} for a deterministic quadrature approximation and \code{method = "MC"} for a stochastic Monte Carlo approximation. See \bold{Details} for more information. 		
}
}

\details{
If \code{d1} is an \code{"ace"} object, \code{family}, \code{formula}, \code{criterion}, \code{prior} and \code{method} will be extracted. If these do not agree with any specified arguments to the function, a warning will be returned

Three different pseudo-Bayesian relative efficiencies are implemented:

\enumerate{
\item
\bold{Pseudo-Bayesian relative D-efficiency} (\code{criterion = "D"})

Here the efficiency is
\deqn{eff^D(d_1,d_2) = 100 \times \exp \left( \left( U^D(d_1) - U^D(d_2) \right)/p \right)}{
eff^D(d_1,d_2) = \exp ( ( U^D(d_1) - U^D(d_2) )/p )}
where
\deqn{U^D(d) = E \left( \log \vert \mathcal{I} (\theta;d) \vert \right)}{
U^D(d) = E ( log | I(\theta;d)|)}
and \eqn{p}{p} is the number of parameters.
\item
\bold{Pseudo-Bayesian A-optimal} (\code{criterion = "A"})

Here the efficiency is
\deqn{eff^A(d_1,d_2) = 100 \times U^A(d_2)/U^A(d_1) }{
eff^A(d_1,d_2) = 100 \times U^A(d_2)/U^A(d_1)}
where
\deqn{U^A(d) = - E \left( \mathrm{tr} \left\{ \mathcal{I} (\theta;d)^{-1} \right\} \right)}{
U^A(d) = - E ( tr (I(\theta;d)^(-1)))}
\item
\bold{Pseudo-Bayesian E-optimal} (\code{criterion = "E"})

Here the efficiency is
\deqn{eff^E(d_1,d_2) = 100 \times U^E(d_2)/U^E(d_1) }{
eff^E(d_1,d_2) = 100 \times U^E(d_2)/U^E(d_1)}
where
\deqn{U^E(d) = E \left( \mathrm{min} \mbox{ } e\left(\mathcal{I} (\theta;d) \right) \right)}{
U^E(d) = E ( min(e(I(\theta;d))))}
and \eqn{e()}{e()} denotes the function that calculates the eigenvalues of its argument.}

The expected utilities can be approximated using Monte Carlo methods (\code{method = "MC"} for all criteria) or using a deterministic quadrature method (\code{method = "quadrature"}, implemented for the D, A and E criteria). The former approach approximates the expected utility via sampling from the prior. The latter approach uses a radial-spherical integration rule (Monahan and Genz, 1997) and \code{B[1]} specifies the number, \eqn{n_r}{n_r}, of radial abscissas and \code{B[2]} specifies the number,  \eqn{n_q}{n_q}, of random rotations. Larger values of  \eqn{n_r}{n_r} will produce more accurate, but also more computationally expensive, approximations. See Gotwalt et al. (2009) for further details.   
}

\value{
The function will return a scalar giving the approximated efficiency.
}

\references{
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
## This example involves finding a Bayesian D-optimal design for a 
## compartmental model with n = 18 runs. There are three parameters. 
## Two parameters have uniform priors and the third has a prior 
## point mass. 

n <- 18
k <- 1
p <- 3
set.seed(1)
start.d <- randomLHS(n = n, k = k) * 24
colnames(start.d) <- c("t")

a1<-c(0.01884, 0.298)
a2<-c(0.09884, 8.298)

prior <- list(support = cbind(rbind(a1, a2), c(21.8, 21.8)))
colnames(prior[[1]]) <- c("theta1", "theta2", "theta3") 

example <- acenlm(formula = ~ theta3 * (exp( - theta1 * t) - exp( - theta2 * t)), 
start.d = start.d, prior = prior, lower = 0, upper = 24, N1 = 2, N2 = 0)

## Compute efficiency of final design compared to starting design.
effnlm(d1=example, d2 = start.d)

## Should get 396.9804

}
