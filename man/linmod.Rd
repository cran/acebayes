\name{linmod}
\alias{utillinmod}
\alias{optdeslinmod}

\title{
Example functions for second-order response surface in two factors
}
\encoding{UTF-8}
\description{
This suite of functions implement the example in the Supplementary Material to Overstall & Woods (2015).}
\usage{
utillinmod(d, B)

optdeslinmod(n, type = "ACE")
}

\arguments{
  \item{d}{
An \code{n} by 2 matrix specifying the design matrix. The elements of \code{d} should be in the interval [-1,1].
}
  \item{B}{
A scalar integer specifying the Monte Carlo sample size.
}
	\item{n}{
An integer specifying the total number of runs in the experiment. This argument can be in the set \eqn{\left\{6,7,8,9\right\}}{{6,7,8,9}}.
}
	\item{type}{
An optional character argument specifying which design to return. Possible values are \code{c("ACE","BoxDraper")}. If \code{"ACE"} (the default) then the design found by the ACE algorithm will be returned. If \code{"BoxDraper"} then the true optimal design, as found by Box & Draper (1971), will be returned.
}
}
\details{
This example considers finding D-optimal designs for a second-order response surface in \code{n=6,7,8,9} runs. Note that the D-optimal design is equivalent to the optimal design under Shannon information gain and a non-informative prior distribution.

The expected utility function in this case is available in closed form, i.e. it does not require approximation. Box & Draper (1971) found optimal designs analytically for the number of runs considered here. Overstall & Woods (2015) use this example to demonstrate the efficacy of the ACE algorithm.

For more details, see the Supplementary Material of Overstall & Woods (2015).
}

\value{
The function \code{utillinmod} will return a vector of length \code{B} where each element is a realisation of a stochastic approximation to the expected utility.

The function \code{optdeslinmod} will return an \code{n} by 2 matrix giving the D-optimal design (specified by the argument \code{type}). If \code{type = "ACE"}, the designs returned by this function are those found using the ACE algorithm and are reported in the Supplementary Material of Overstall & Woods (2015), and if \code{type = "BoxDraper"}, the designs returned are the exact D-optimal designs.
}

\references{
Box, M. & Draper, N. (1971). Factorial designs, the \eqn{|F^T F|} criterion and some related matters. 
\emph{Techometrics}, \bold{13}, 731-742.

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
## Set seed for reproducibility

d<-cbind(rep(c(-1,0,1),each=3),rep(c(-1,0,1),3))
## Create a 9-run design which is the true D-optimal design

utillinmod(d=d, B=5)
## Calculate the approximation to the true expected D-optimal utility 
## function for a sample of size 5. 

#[1]  7.926878  8.736976  7.717704 10.148613  8.882840

d0<-optdeslinmod(9)
## Optimal D-optimal design found using ACE

X<-cbind(1,d,d^2,d[,1]*d[,2])
X0<-cbind(1,d0,d0^2,d0[,1]*d0[,2])
## Calculate model matrices under both designs

detX<-determinant(t(X)\%*\%X)$modulus[1]
detX0<-determinant(t(X0)\%*\%X0)$modulus[1]
## Calculate true expected D-optimal utility function for both designs

100*exp(0.2*(detX0-detX))
## Calculate D-efficiency of ACE design.

# 99.93107
}
