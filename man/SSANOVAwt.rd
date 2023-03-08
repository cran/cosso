\name{SSANOVAwt}
\alias{SSANOVAwt}
\title{
Compute adaptive weights by fitting a SS-ANOVA model
}


\description{
A preliminary estimate \eqn{\tilde{\eta}} is first obtained by fitting a smoothing spline ANOVA model,
and then use the inverse \eqn{L_2}-norm, \eqn{||\tilde{\eta}_j||^{-\gamma}}, as the initial weight for the \eqn{j}-th functional component.
}

\usage{ SSANOVAwt(x,y,tau,family=c("Gaussian","Binomial","Cox","Quantile"),mscale=rep(1,ncol(x)),
               gamma=1,scale=FALSE,nbasis,basis.id,cpus) }

\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension.
         The range of input variables is scaled to [0,1] for continuous variables.}
\item{y}{response vector. Quantitative for \code{family="Gaussian"} or \code{family="Quantile"}.
         For \code{family="Binomial"} should be a vector with two levels.
         For \code{family="Cox"}, y should be a two-column matrix (data frame) with columns named 'time' and 'status'}
\item{tau}{the quantile to be estimated, a number strictly between 0 and 1. Argument required when \code{family="Quantile"}.}
\item{family}{response type. Abbreviations are allowed.}
\item{mscale}{scale parameter for the Gram matrix associated with each function component. Default is \code{rep(1,ncol(x))}}
\item{gamma}{power of inverse \eqn{L_2}-norm. Default is \code{1}.}
\item{scale}{if \code{TRUE}, continuous predictors will be rescaled to [0,1] interval. Default is \code{FALSE}.}
\item{nbasis}{number of "knots" to be selected. Ignored when \code{basis.id} is provided.}
\item{basis.id}{index designating selected "knots". Argument is not valid if \code{family="Quantile"}.}
\item{cpus}{number of available processor units. Default is \code{1}. If \code{cpus}>=2, parallelize task using "parallel" package. Recommended when either sample size or number of covariates is large.
            Argument is not valid if \code{family="Gaussian"} or \code{family="Binomial"}.}
}

\details{
The initial mean function is estimated via a smooothing spline objective function. In the SS-ANOVA model framework,
the regression function is assumed to have an additive form
\deqn{\eta(x)=b+\sum_{j=1}^p\eta_j(x^{(j)}),}
where \eqn{b} denotes intercept and \eqn{\eta_j} denotes the main effect of the \eqn{j}-th covariate.

For \code{"Gaussian"} response, the mean regression function is estimated by minimizing the objective function:
\deqn{\sum_i(y_i-\eta(x_i))^2/nobs+\lambda_0\sum_{j=1}^p\alpha_j||\eta_j||^2.}
where RSS is residual sum of squares.

For \code{"Binomial"} response, the regression function is estimated by minimizing the objective function:
\deqn{-log-likelihood/nobs+\lambda_0\sum_{j=1}^p\alpha_j||\eta_j||^2}

For \code{"Quantile"} regression model, the quantile function, is estimated by minimizing the objective function:
\deqn{\sum_i\rho(y_i-\eta(x_i))/nobs+\lambda_0\sum_{j=1}^p\alpha_j||\eta_j||^2.}

For \code{"Cox"} regression model, the log-hazard function, is estimated by minimizing the objective function:
\deqn{-log-Partial Likelihood/nobs+\lambda_0\sum_{j=1}^p\alpha_j||\eta_j||^2.}

The smoothing parameter \eqn{\lambda_0} is tuned by 5-fold Cross-Validation, if \code{family="Gaussian"}, \code{"Binomial"} or \code{"Quantile"},
and Approximate Cross-Validation, if \code{family="Cox"}. But the smoothing parameters \eqn{\alpha_j} are given in the argument \code{mscale}.

The adaptive weights are then fiven by \eqn{||\tilde{\eta}_j||^{-\gamma}}.
}

\value{
\item{wt}{The adaptive weights.}
}

\references{
Storlie, C. B., Bondell, H. D., Reich, B. J. and Zhang, H. H. (2011) "Surface Estimation, Variable Selection, and the Nonparametric Oracle Property", Statistica Sinica, \bold{21}, 679--705.

}

\author{
Hao Helen Zhang and Chen-Yen Lin}


\examples{
## Adaptive COSSO Model
## Binomial
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*7,0,1),nc=7))
trueProb=1/(1+exp(-x[,1]-sin(2*pi*x[,2])-5*(x[,4]-0.4)^2))
y=rbinom(200,1,trueProb)

Binomial.wt=SSANOVAwt(x,y,family="Bin")
ada.B.Obj=cosso(x,y,wt=Binomial.wt,family="Bin")

\dontrun{
## Gaussian
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*7,0,1),nc=7))
y=x[,1]+sin(2*pi*x[,2])+5*(x[,4]-0.4)^2+rnorm(200,0,1)
Gaussian.wt=SSANOVAwt(designx,response,family="Gau")
ada.G.Obj=cosso(x,y,wt=Gaussian.wt,family="Gaussian")
}

}
