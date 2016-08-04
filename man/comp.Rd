\name{comp}
\alias{utilcomp18bad}
\alias{optdescomp18bad}
\alias{utilcomp15bad}
\alias{optdescomp15bad}
\alias{optdescomp15sig}
\alias{utilcomp15sig}
\alias{utilcomp15sigDRS}
\alias{optdescomp15sigDRS}

\title{
Compartmental model.
}
\encoding{UTF-8}
\description{
This suite of functions implement the examples in \emph{Section 3.2: Compartmental model} in Overstall & Woods (2015) and the accompanying Supplementary Material.
}
\usage{
utilcomp18bad(d, B)

optdescomp18bad(type = "ACE")

utilcomp15bad(d, B)

optdescomp15bad()

utilcomp15sig(d, B)

optdescomp15sig()

utilcomp15sigDRS(d, B)

optdescomp15sigDRS()

}

\arguments{
  \item{d}{
An \code{n} by 1 matrix specifying the design matrix controlling the sampling times. Here \code{n} can be in the set \eqn{\left\{2,15,18\right\}}{{2,15,18}} depending on the application (see \bold{Details} and \bold{Value}). If \code{n} is 15 or 18, each element should be scaled to be in the interval [-1, 1]; otherwise the elements should be postive.
}
  \item{B}{
A scalar integer specifying the Monte Carlo sample size.
}
  \item{type}{
An optional character argument specifying which design to return. Possible values are \code{c("ACE","Gotwalt","Atkinson")}. If \code{"ACE"} (the default) then the design found by the ACE algorithm will be returned. If \code{"Gotwalt"} then the design published in Gotwalt et al (2009) is returned. If \code{"Atkinson"} then the design found by Atkinson et al (1993) is returned.}

}
\details{
Compartmental models are used in Pharmokinetics to study how materials 
ow through an organism. A drug is administered to an individual or animal and then the 
amount present at a certain body location is measured at a set of n pre-determined sampling 
times (in hours). There is one design variable: sampling time. Therefore the design matrix is 
an \code{n} by 1 matrix with elements controlling the \code{n} sampling times.

In Overstall & Woods (2015), two different compartmental model examples are considered. The first (in the Supplementary Material) 
comes from Atkinson et al (1993) and Gotwalt et al (2009) where there are \code{n = 18} sampling times and interest 
lies in finding a Bayesian D-optimal design. The functions whose name includes \code{"comp18"} refer to this example.

The second example (in Section 3.2) comes from Ryan et al (2014), where there are \code{n=15} sampling times and 
the ultimate interest lies in finding an optimal design under the Shannon information gain utility. Also considered is the
Bayesian D-optimal design. The functions whose name includes \code{"comp15"} refer to this example. Note that Ryan et al (2014) used a dimension reduction scheme (DRS) to find optimal designs. The function whose name is suffixed by \code{"DRS"} refer to this situation.
}

\value{
For the example in the Supplementary Material;
\itemize{
\item{
The function \code{utilcomp18bad} will return a vector of length \code{B} where each element is the value of the Bayesian D-optimal utility function (i.e. the log determinant of the Fisher information matrix) evaluated at a sample of size \code{B} generated from the prior distribution of the model parameters.}
\item{
The function \code{optdescomp18bad} will return an 18 by 1 matrix giving the optimal design (specified by the argument \code{type}). The elements will be scaled to be in the interval [-1, 1], i.e. a -1 corresponds to a sampling times of 0 hours, and 1 corresponds to a time of 24 hours.}}
For the example in Section 3.2;
\itemize{
\item{
The function \code{utilcomp15bad} will return a vector of length \code{B} where each element is the value of the Bayesian D-optimal utility function (i.e. the log determinant of the Fisher information matrix) evaluated at a sample of size \code{B} generated from the prior distribution of the model parameters.}
\item{
The function \code{optdescomp15bad} will return an 15 by 1 matrix giving the optimal design (found using ACE) under Bayesian D-optimality. The elements will be scaled to be in the interval [-1, 1], i.e. a -1 corresponds to a sampling times of 0 hours, and 1 corresponds to a time of 24 hours.}
\item{
The function \code{utilcomp15sig} will return a vector of length \code{B} where each element is the value of the SIG utility function evaluated at a sample of size \code{B} generated from the joint distribution of model parameters and unobserved responses.}
\item{
The function \code{optdescomp15sig} will return an 18 by 1 matrix giving the optimal design (found using ACE) under the SIG utility. The elements will be scaled to be in the interval [-1, 1], i.e. a -1 corresponds to a sampling times of 0 hours, and 1 corresponds to a time of 24 hours.}
\item{
The function \code{utilcomp15sigDRS} will return a vector of length \code{B} where each element is the value of the SIG utility function (where a DRS has been used) evaluated at a sample of size \code{B} generated from the joint distribution of model parameters and unobserved responses. Here the Beta DRS (see Overstall & Woods, 2015) has been used so \code{d} should be a 2 by 1 matrix containing the positive beta parameters.}
\item{
The function \code{optdescomp15sigDRS} will return a 2 by 1 matrix giving the optimal design (found using ACE) under the SIG utility, where a DRS has been used. The elements correspond to the parameters of a beta distribution.}
}}

\references{
Atkinson, A., Chaloner, K., Herzberg, A., & Juritz, J. (1993). Experimental Designs for Properties
of a Compartmental Model. \emph{Biometrics}, \bold{49}, 325-337.

Gotwalt, C., Jones, B. & Steinberg, D. (2009). Fast Computation of Designs Robust to Parameter
Uncertainty for Nonlinear Settings. \emph{Technometrics}, \bold{51}, 88-95.

Overstall, A.M. & Woods, D.C. (2015). Bayesian Design of Experiments using
Approximate Coordinate Exchange. \emph{arXiv:1501.00264}.

Ryan, E., Drovandi, C., Thompson, M., Pettitt, A. (2014). Towards Bayesian experimental design
for nonlinear models that require a large number of sampling times. \emph{Computational
Statistics and Data Analysis}, \bold{70}, 45-60.
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

d<-optimumLHS(n=18,k=1)*2-1
## Generate an 18-run design.

u<-utilcomp18bad(d=d, B=20000)
## Calculate the D-optimal utility function for a 
## sample of size 20000. 

u[1:5]
## Look at first 5 elements.
#[1] 14.283473 10.525390  4.126233  7.061498 12.793245

d0<-optdescomp18bad()
u0<-utilcomp18bad(d=d0, B=20000)
## Optimal design found by ACE and calculate the D-optimal 
## utility function for a sample of size 20000.

u0[1:5]
## Look at first 5 elements.
#[1] 15.04721 15.37141 16.84287 14.06750 14.01523

mean(u)
mean(u0)
## Caluclate expected Bayesian D-optimal utility.

d<-matrix(runif(2),ncol=1)
## Generate two beta parameters.

u<-utilcomp15sigDRS(d=d, B=5)
u
## Calculate the SIG utility function for a 
## sample of size 5.
#[1] 17.652044  4.878998 19.919559 22.017760  5.600473
}
