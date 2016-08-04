\name{beetle}
\alias{utilbeetle}
\alias{optdesbeetle}

\title{
Example functions for binomial regression under model uncertainty
}
\encoding{UTF-8}
\description{
This suite of functions implement the example in \emph{Section 3.4: Binomial regression under model uncertainty} in Overstall & Woods (2015). 
}
\usage{
utilbeetle(d, B)

optdesbeetle(n)
}

\arguments{
  \item{d}{
An \code{n} by 1 matrix specifying the set of \code{n} new doses. The elements of \code{d} should be scaled to be in the interval [-1, 1].
}
  \item{B}{
A scalar integer specifying the Monte Carlo sample size.
}
  \item{n}{
A scalar integer specifying the number of new doses. \code{n} should be in the set \eqn{\left\{1,\dots,10\right\}}{{1,...,10}}.
}
}
\details{
This example considers generating an optimal design for a follow-up experiment. The original design and data (Bliss, 1935) involves administering different doses of poison to \eqn{N = 8} groups of beetles. The number of beetles that die in each group are recorded. Six different models are considered formed from the combination of three link functions and two linear predictors. Interest lies in the quantity known as lethal dose 50 (LD50) which is the dose required to kill 50\% of the beetles and is a function of the model parameters for a given model. Consider finding an optimal design for estimating LD50 under the negative squared error loss (NSEL) function for \code{n} new doses of poison. The prior distribution is equivalent to the posterior distribution arising from the original data and includes model uncertainty.

For more details, Overstall & Woods (2015).
}

\value{
The function \code{utilbeetle} will return a vector of length \code{B} where each element is the value of the utility function for a sample generated from the joint distribution of the model parameters, model and unobserved responses.

The function \code{optdebeetle} will return an \code{n} by 1 matrix giving the optimal design under the NSEL utility function (found using ACE) for estimating the LD50. The elements will be scaled to be in the interval [-1, 1], where -1 corresponds to a dose of 1.6907, 0 to a dose of 1.7873 and 1 to a dose of 1.8839. The designs given by this function are those reported in Overstall & Woods (2015) for \code{n} = 1, \ldots, 10.
}

\references{
Bliss, C. (1935). The calculation of the dosage-mortality curve. \emph{Annals of Applied Biology},
\bold{22}, 134-167.

O'Hagan, A, & Forster, J.J. (2004). \emph{Kendall's Advanced Theory of Statistics, Volume 2B: Bayesian
Inference. 2nd edition}. John Wiley & Sons.

Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.
}

\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk} & David C. Woods.
}

\seealso{
\code{\link{ace}}.
}
\examples{
set.seed(1)
## Set seed for repoducibility

d<-optimumLHS(n=10,k=1)*2-1
## Generate a design of 10 doses with elements in [-1, 1]

utilbeetle(d=d, B=5)
## Calculate the utility function for a sample of size 5.

#-4.085626e-06 -1.209499e-06 -1.214244e-05 -6.243750e-06 -3.083291e-05

d0<-optdesbeetle(10)
d0
## Print out optimal design from Overstall & Woods (2015) for 10 doses 

0.5*(d0+1)*(1.8839-1.6907)+1.6907
## On original scale.
#          [,1]
# [1,] 1.769957
# [2,] 1.769520
# [3,] 1.768641
# [4,] 1.777851
# [5,] 1.768641
# [6,] 1.769520
# [7,] 1.777851
# [8,] 1.765997
# [9,] 1.768641
#[10,] 1.768641

}
