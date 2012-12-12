\name{SSANOVAwt.cox}
\alias{SSANOVAwt.cox}
\title{
Compute adaptive weights for Cox proportional hazards regression} 

\description{
A preliminary estimate of log-relative risk \eqn{\eta} is first obtained by minimizing a smoothing
spline type of objective function, 
and then the weight for the j-th function component is computed as \eqn{||P_j \tilde{\eta}||^{-1}}. 
Here \eqn{P_j} denotes the projection operator to the subspace. 
By default, we use the reciprocal \eqn{L_2} norm as the initial weight.
}

\usage{ SSANOVAwt.cox(x,time,status,mscale=rep(1,ncol(x))) }


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension.
         The range of input variables is scaled to [0,1].}
\item{time}{observed event time.}
\item{status}{status indicator. 1: dead; 0: alive.}
\item{mscale}{scale parameter for the Gram matrix associated with each function component. Default is \code{rep(1,ncol(x))}}
}

\details{
The initial log-relative risk is estimated via the penalized partial likelihood:
\deqn{-partial~likelihood/nobs+\lambda_0\sum_{j=1}^p\alpha_j||P_j\eta||^2.}

The smoothing parameters \eqn{\lambda_0} will be tuned by approximate cross-validation (Leng and Zhang 2006),
but \eqn{\alpha_j}'s are given in the argumnent \code{mscale}.
}


\value{
\item{adawt}{The adaptive weights}
}


\author{
Hao Helen Zhang and Chen-Yen Lin}


\examples{
data(veteran)
adawt <- SSANOVAwt.cox(x=veteran[,c(3,5:8)],time=veteran[,1],status=veteran[,2])
print(adawt)
}
