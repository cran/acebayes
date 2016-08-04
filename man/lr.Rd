\name{lr}

\alias{utillrbad}
\alias{optdeslrbad}
\alias{utillrsig}
\alias{inideslrsig}
\alias{optdeslrsig}

\alias{utillrbaa}
\alias{optdeslrbaa}
\alias{utillrnsel}
\alias{inideslrnsel}
\alias{optdeslrnsel}

\alias{utilhlrbad}
\alias{optdeshlrbad}
\alias{utilhlrsig}
\alias{inideshlrsig}
\alias{optdeshlrsig}

\alias{utilhlrbaa}
\alias{optdeshlrbaa}
\alias{utilhlrnsel}
\alias{inideshlrnsel}
\alias{optdeshlrnsel}

\title{
Logistic regression in four factors.
}
\encoding{UTF-8}
\description{
This suite of functions implement the examples in \emph{Section 3.3: Logistic regression in four factors} in Overstall & Woods (2015).
}
\usage{
utillrbad(d, B)
optdeslrbad(n, type = "ACE")
utillrsig(d, B)
inideslrsig(n, rep)
optdeslrsig(n)

utilhlrbad(d, B)
optdeshlrbad(n)
utilhlrsig(d, B)
inideshlrsig(n, rep)
optdeshlrsig(n)

utillrbaa(d, B)
optdeslrbaa(n)
utillrnsel(d, B)
inideslrnsel(n, rep)
optdeslrnsel(n)

optdeshlrbaa(n)
utilhlrbaa(d, B)
utilhlrnsel(d, B)
inideshlrnsel(n, rep)
optdeshlrnsel(n)
}

\arguments{
  \item{d}{
An \code{n} by 4 matrix specifying the design matrix. Each element should be in the interval [-1, 1]. The value of \code{n} depends on the application (see \bold{Value}).
}
  \item{B}{
A scalar integer specifying the Monte Carlo sample size.
}
  \item{n}{
A scalar integer specifying the number of runs whose value depends on the application (see \bold{Value}).
}
  \item{rep}{
A scalar integer in the set \eqn{\left\{1,\dots,20\right\}}{{1,...,20}} specifying the initial design. 
}
  \item{type}{
An optional character argument specifying which design to return. Possible values are \code{c("ACE","Gotwalt1","Gotwalt2","Woods")}. If \code{"ACE"} (the default) then the design found by the ACE algorithm will be returned. If \code{"Gotwalt1"} then the design published in Gotwalt et al (2009) is returned (only available for \code{n = 16}). If \code{"Gotwalt2"} then the design found by Overstall & Woods (2015) using the Gotwalt et al (2009) method is returned (only available for \code{n = 16, 48}). If \code{"Woods"} then the design found by Woods et al (2006) are returned (only available for \code{n = 16, 48}).}
}
\details{
This example involves a first-order logistic regression model in four factors (and \code{n} runs). Woods et al (2006) and Gotwalt et al (2009) considered generating Bayesian D-optimal designs for \code{n = 16} and \code{n = 48}. Overstall & Woods (2015) extended this example by considering Bayesian A-optimal, Shannon information gain (SIG) and negative sqaured error loss utility functions, a range of number of runs from 6 to 48, and "random effects" to form a hierarchical logistic regression model.

See Overstall & Woods (2015) for more details on this example.
}

\value{
A function whose name includes \code{"lr"} refers to standard logistic regression, whereas \code{"hlr"} refers to hierarchical logistic regression. Under standard logistic regression the possible values for the argument \code{n} can be any even integer between 6 and 48. For hierarchical logistic regression, \code{n} can be any integer divisible by 6 between 12 and 48. The function name also indicates the utility function:
\itemize{
\item{\code{"bad"}}{
Bayesian D-optimal}
\item{\code{"baa"}}{
Bayesian A-optimal}
\item{\code{"sig"}}{
Shannon information gain}
\item{\code{"nsel"}}{
Negative squared error loss}}
The functions prefixed by \code{"util"} will return a vector of length \code{B} where each element is the utility function evaluated at a sample generated from the prior distribution of model parameters (for Bayesian D- and A-optimality) or the joint distribution of model parameters and unobserved responses (for SIG and NSEL).

The functions prefixed by \code{"optdes"} will return an \code{n} by 4 matrix giving the optimal design found by ACE. The designs given by this function are those reported on in Overstall & Woods (2015). The function \code{optdeslrbad} will also return designs (for \code{n = 16, 48}) found by Woods et al (2006) and Gotwalt et al (2009) by specifying the \code{type} argument appropriately.

The functions prefixed by \code{"inides"} will return an \code{n} by 4 matrix giving an initial design for ACE to find the optimal designs under the SIG and NSEL utility functions. These are 20 designs found using ACE under approximations to the Bayesian A- and D-optimal utility functions, respectively. The argument \code{rep} specifies which of these 20 designs to use.
}

\references{
Gotwalt, C., Jones, B. & Steinberg, D. (2009). Fast Computation of Designs Robust to Parameter
Uncertainty for Nonlinear Settings. \emph{Technometrics}, \bold{51}, 88-95.

Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.

Woods, D.C., Lewis, S., Eccleston, J., Russell, K. (2006). Designs for Generalized Linear Models
With Several Variables and Model Uncertainty. \emph{Technometrics}, \bold{48}, 284-292.
}

\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk} & David C. Woods.
}

\seealso{
\code{\link{ace}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility

d<-optimumLHS(n=16,k=4)*2-1
## Generate an 16-run design.

u<-utillrbad(d=d, B=20000)
## Calculate the D-optimal utility function for a 
## sample of size 20000. 

u[1:5]
## Look at first 5 elements.
#[1] -11.426457  -5.925467  -9.238526  -9.392765  -7.689553

d0<-optdeslrbad(16)
u0<-utillrbad(d=d0, B=20000)
## Optimal design found by ACE and calculate the D-optimal 
## utility function for a sample of size 20000.

u0[1:5]
## Look at first 5 elements.
#[1] -4.644116 -2.411431 -4.999891 -2.906558 -2.282687

mean(u)
mean(u0)
## Calculate expected Bayesian D-optimal utility.
#[1] -9.349215
#[1] -2.992012
}
