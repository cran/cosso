\name{compwt}
\alias{compwt}
\title{
Computing the weights used in the penalty term of the adaptive COSSO.} 

\description{
A preliminary estimate of \eqn{f} is first obtained, 
and then the weight for the jth function component is computed as \eqn{||P_j f||^{-\gamma}}. 
Here \eqn{P_j} denotes the projection operator to the subspace. By default, we use the smoothing spline solution 
as the initial estimator and set \eqn{\gamma} as 1.
}

\usage{ compwt(x,y,mscale=rep(1,ncol(x)),gampow=1) }


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension.
         The range of input variables is scaled to [0,1].}
\item{y}{response vector}
\item{mscale}{Scale parameter for the Gram matrix associated with each function component. Default is \code{rep(1,ncol(x))}}
\item{gampow}{Power of the \eqn{L_2} norm. Default is \code{1}}
}


\value{
\item{adwt}{The adaptive weights used in the adaptive COSSO}
}

\references{Storlie, C, Bondell, H., Reich, B. and Zhang, H.H. "Surface estimation, variable selection, and the nonparametric oracle property", 
Statistica Sinica, 2011.}

\author{
Hao Helen Zhang \email{hzhang@stat.ncsu.edu} }


\examples{
data(ozone)
set.seed(2010)
# Randomly select one-third of the data as training data.
train_id <- sample(1:nrow(ozone),round(nrow(ozone)/3))
wt <- compwt(ozone[train_id,-1],ozone[train_id,1],rep(1,ncol(ozone)-1),1)
ozone_cosso <- cosso(x=ozone[train_id,-1],y=ozone[train_id,1],wt=wt,type="BIC")
summary.cosso(ozone_cosso)
}
