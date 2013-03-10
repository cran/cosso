\name{cosso} 
\alias{cosso}
\title{Fit a generalized nonparametric model with cosso penalty}

\description{
A comprehensive method for fitting various type of regularized nonparametric regression models using cosso penalty.
Fits mean, logistic, Cox and quantile regression.}

\usage{
cosso(x,y,tau,family=c("Gaussian","Binomial","Cox","Quantile"),wt=rep(1,ncol(x)),
      scale=FALSE,nbasis,basis.id,cpus)
}


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension. 
         The range of input variables is scaled to [0,1] for continuous variables. Variables with less than 7 unique values will be considered as discrete variable.}
\item{y}{response vector. Quantitative for \code{family="Gaussian"} or \code{family="Quantile"}. 
         For \code{family="Binomial"} should be a vector with two levels.
         For \code{family="Cox"}, y should be a two-column matrix (or data frame) with columns named 'time' and 'status'}
\item{tau}{the quantile to be estimated, a number strictly between 0 and 1. Arguement required when \code{family="Quantile"}.}
\item{family}{response type. Abbreviations are allowed.}
\item{wt}{weights for predictors. Default is \code{rep(1,ncol(x))}}
\item{scale}{if \code{TRUE}, continuous predictors will be rescaled to [0,1] interval. Dafault is \code{FALSE}.}
\item{nbasis}{number of "knots" to be selected. Ignored when \code{basis.id} is provided.}
\item{basis.id}{index designating selected "knots". Arguement is not valid for \code{family="Quantile"}.}
\item{cpus}{number of available processor units. Default is \code{1}. If \code{cpus}>=2, parallelize task using "parallel" package. Recommended when either sample size or number of covariates is large.
            Arguement is only valid for \code{family="Cox"} or \code{family="Quantile"}.}
}

\details{
In the SS-ANOVA model framework, the regression function is assumed to have an additive form
\deqn{\eta(x)=b+\sum_{j=1}^p\eta_j(x^{(j)}),}
where \eqn{b} denotes intercept and \eqn{\eta_j} denotes the main effect of the \eqn{j}-th covariate.

For \code{"Gaussian"} response, the mean function is estimated by minimizing the objective function:
\deqn{\sum_i(y_i-\eta(x_i))^2/nobs+\lambda_0\sum_{j=1}^p\theta^{-1}_jw_j^2||\eta_j||^2, s.t.~\sum_{j=1}^p\theta_j\leq M.}

For \code{"Binomial"} response, the log-odd function is estimated by minimizing the objective function:
\deqn{-log-likelihood/nobs+\lambda_0\sum_{j=1}^p\theta^{-1}_jw_j^2||\eta_j||^2, s.t.~\sum_{j=1}^p\theta_j\leq M.}

For \code{"Quantile"} regression model, the quantile function, is estimated by minimizing the objective function:
\deqn{\sum_i\rho_{\tau}(y_i-\eta(x_i))/nobs+\lambda_0\sum_{j=1}^p\theta^{-1}_jw_j^2||\eta_j||^2, s.t.~\sum_{j=1}^p\theta_j\leq M.}

For \code{"Cox"} regression model, the log-relative hazard function is estimated by minimizing the objective function:
\deqn{-log-Partial Likelihood/nobs+\lambda_0\sum_{j=1}^p\theta^{-1}_jw_j^2||\eta_j||^2, s.t.~\sum_{j=1}^p\theta_j\leq M.}
For identifiability sake, the intercept term in Cox model is absorbed into basline hazard, or equivalently set \eqn{b=0}.


For large data sets, we can reduce the computational load of the optimization problem by
selecting a subset of the observations of size \emph{nbais} as "knots", which reduces the dimension of the
kernel matrices from \emph{nobs} to \emph{nbasis}.
Unless specified via \code{basis.id} or \code{nbasis}, the default number of "knots" is \eqn{max}(40,12*nobs^(2/9)) for \code{"Gaussian"} and \code{"Binomial"} and
\eqn{max}(35,11 * nobs^(2/9)) for \code{"Cox"}.


The weights can be specified based on either user's own discretion or adaptively computed from initial 
function estimates. See Storlie et al. (2011) for more discussions. One possible choice is to specify the weights
as the inverse \eqn{L_2} norm of initial function estimator, see \code{\link{SSANOVAwt}}.
}

\value{
An object with S3 class "cosso".
\item{y}{the response vector.}
\item{x}{the input matrix.}
\item{Kmat}{a three-dimensional array containing kernel matrices for each input variables.}
\item{wt}{weights for predictors.}
\item{family}{type of regression model.}
\item{basis.id}{indices of observations used as "knots".}
\item{cpus}{number of cpu units used. Will be returned if \code{family="Cox"} or \code{family="Quantile"}.}
\item{tau}{the quantile to be estimated. Will be returned if \code{family="Quantile"}.}
\item{tune}{a list containing prelminary tuning result and L2-norm.}
}

\references{
Lin, Y. and Zhang, H. H. (2006) "Component Selection and Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.

Leng, C. and Zhang, H. H. (2006) "Model selection in nonparametric hazard regression", Nonparametric Statistics, \bold{18}, 417--429.

Zhang, H. H. and Lin, Y. (2006) "Component Selection and Smoothing for Noparametric Regression in Exponential Families", Statistica Sinica, \bold{16}, 1021--1041.

Storlie, C. B., Bondell, H. D., Reich, B. J. and Zhang, H. H. (2011) "Surface Estimation, Variable Selection, and the Nonparametric Oracle Property", Statistica Sinica, \bold{21}, 679--705.

}

\author{
Hao Helen Zhang and Chen-Yen Lin}

\seealso{ \code{\link{plot.cosso}}, \code{\link{predict.cosso}}, \code{\link{tune.cosso}}
}

\examples{
## Gaussian
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*9,0,1),nc=9))
y=x[,1]+sin(2*pi*x[,2])+5*(x[,4]-0.4)^2+rnorm(200,0,1)
G.Obj=cosso(x,y,family="Gaussian")
plot.cosso(G.Obj,plottype="Path")

\dontrun{
## Use all observations as knots
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*9,0,1),nc=9))
y=x[,1]+sin(2*pi*x[,2])+5*(x[,4]-0.4)^2+rnorm(200,0,1)
G.Obj=cosso(x,y,family="Gaussian",nbasis=200)
## Clean up
rm(list=ls())

## Binomial
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*9,0,1),nc=9))
trueProb=1/(1+exp(-x[,1]-sin(2*pi*x[,2])-5*(x[,4]-0.4)^2))
y=rbinom(200,1,trueProb)

B.Obj=cosso(x,y,family="Bin") 
## Clean up
rm(list=ls())

## Cox
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*9,0,1),nc=9))
hazard=x[,1]+sin(2*pi*x[,2])+5*(x[,4]-0.4)^2
surTime=rexp(200,exp(hazard))
cenTime=rexp(200,exp(-hazard)*runif(1,4,6))
y=cbind(time=apply(cbind(surTime,cenTime),1,min),status=1*(surTime<cenTime))
C.obj=cosso(x,y,family="Cox",cpus=1)

## Try parallel computing
C.obj=cosso(x,y,family="Cox",cpus=4)

## Clean up
rm(list=ls())

## Quantile
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*7,0,1),nc=7))
y=x[,1]+sin(2*pi*x[,2])+5*(x[,4]-0.4)^2+rt(200,3)
Q.obj=cosso(x,y,0.3,family="Quan",cpus=1)

## Try parallel computing
Q.obj=cosso(x,y,0.3,family="Quan",cpus=4)
}

}
