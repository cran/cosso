\name{tune.cosso}
\alias{tune.cosso}
\title{Tuning procedure for cosso}

\description{
Compute K-fold cross-validated score and plot cross-validated score against a grid values of smooth parameter M. 
}

\usage{ 
tune.cosso(object,folds=5,plot.it=TRUE)   
}


\arguments{
\item{object}{a cosso object.}
\item{folds}{number of folds for corss-validation. Default is \code{5}. It is not recommended to use \code{folds} less than 4.}
\item{plot.it}{if \code{TRUE}, plot the cross-validated score against a sequence values of M.}
}


\value{
\item{OptM}{the optimal smoothing parameter for M.}
\item{M}{used tuning grid points.}
\item{cvm}{the mean cross-validated error/minus log-likelihood.}
\item{cvsd}{estimate of standard error of \code{cvm}.}
}


\author{
Hao Helen Zhang and Chen-Yen Lin}

\seealso{\code{\link{cosso}}, \code{\link{predict.cosso}}
}



\examples{ 
## Binomial
set.seed(20130310)
x=cbind(rbinom(150,1,.7),matrix(runif(150*5,0,1),nc=5))
trueProb=1/(1+exp(-x[,1]-sin(2*pi*x[,2])-5*(x[,4]-0.4)^2))
y=rbinom(150,1,trueProb)

B.Obj=cosso(x,y,family="Bin",nbasis=30) 
tune.cosso(B.Obj,4,TRUE)
}
