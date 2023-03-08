\name{predict.cosso} 
\alias{predict.cosso}

\title{
Make predictions or extract coefficients from a cosso model}

\description{
Make prediction for future observations or extract the model parameters at a particular smoothing parameter.
}

\usage{ 
\method{predict}{cosso}(object,xnew,M,type=c("fit","coefficients","nonzero"),eps=1e-7,...) 

}


\arguments{
\item{object}{a cosso object.}
\item{xnew}{matrix of new values for \eqn{x} at which predictions are to be made. Object must be a matrix and have the same dimension as the training design points. Continuous variable will also have to be scaled to [0,1] interval.}
\item{M}{a smoothing parameter value. M should be taken between 0 and p. If not provided, a cross-validation procedure will be carried out to select an appropriate value.}
\item{type}{if \code{type="fit"} (default), fitted values will be returned. If \code{type="coefficients"}, model coefficients will be returned. Abbreviations are allowed.}
\item{eps}{an effective zero, default is \code{1e-7}}
\item{...}{additional arguments for predict function.}
}


\value{
The object returned depends on type.

When \code{type="fit"}, predicted \deqn{eta} function value will be given at the new design points \code{xnew}.

When \code{type="coefficients"}, three sets of coefficients will be returned.
\item{Intercept}{the estimated intercept. If \code{family="Cox"}, the intercept is zero.}
\item{coefs}{the estimated coefficients for kernel representers.}
\item{theta}{the estimated scale parameters for each functional component.}

When \code{type="nonzero"}, a list of the indices of the nonconstant functional components will be returned.
}


\author{
Hao Helen Zhang and Chen-Yen Lin}

\seealso{\code{\link{plot.cosso}}
}

\examples{ 
## Gaussian
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*7,0,1),nc=7))
y=x[,1]+sin(2*pi*x[,2])+5*(x[,4]-0.4)^2+rnorm(200,0,1)
G.Obj=cosso(x,y,family="Gaussian")
predict.cosso(G.Obj,M=2,type="nonzero")
predict.cosso(G.Obj,xnew=x[1:3,],M=2,type="fit")
## Clean up
rm(list=ls())

\dontrun{
## Binomial
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*9,0,1),nc=9))
trueProb=1/(1+exp(-x[,1]-sin(2*pi*x[,2])-5*(x[,4]-0.4)^2))
y=rbinom(200,1,trueProb)

B.Obj=cosso(x,y,family="Bin") 
f.hat=predict.cosso(B.Obj,xnew=x,M=2,type="fit")
prob.hat=1/(1+exp(-f.hat))
## Clean up
rm(list=ls())
}
}
