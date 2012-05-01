\name{SSANOVAwt}
\alias{SSANOVAwt}
\title{
Computing the weights used in the penalty term of the adaptive COSSO.} 

\description{
A preliminary estimate of \eqn{f} is first obtained, 
and then the weight for the jth function component is computed as \eqn{||P_j f||^{-\gamma}}. 
Here \eqn{P_j} denotes the projection operator to the subspace. By default, we use the smoothing spline solution 
as the initial estimator and set \eqn{\gamma} as 1.
}

\usage{ SSANOVAwt(x,y,mscale=rep(1,ncol(x)),gampow=1) }


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension.
         The range of input variables is scaled to [0,1].}
\item{y}{response vector}
\item{mscale}{scale parameter for the Gram matrix associated with each function component. Default is \code{rep(1,ncol(x))}}
\item{gampow}{power of the \eqn{L_2} norm. Default is \code{1}}
}


\value{
\item{adwt}{The adaptive weights used in the adaptive COSSO}
}

\references{
Storlie, C. B., Bondell, H. D., Reich, B. J. and Zhang, H. H. (2011) "Surface Estimation, Variable Selection, and the Nonparametric Oracle Property", Statistica Sinica, \bold{21}, 679--705.

}

\author{
Hao Helen Zhang \email{hzhang@stat.ncsu.edu} }


\examples{
data(ozone)
wt <- SSANOVAwt(ozone[,-1],ozone[,1])
print(wt)
}
