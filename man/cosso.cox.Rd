\name{cosso.cox} 
\alias{cosso.cox}
\title{
Estimate the log-relative risk using Smoothing Splines with COSSO penalty.
}

\description{
Estimate the log-relative risk function in a Cox proportional hazards model
by solving a penalized partial likelihood function.
}

\usage{
cosso.cox(x,time,status,wt=rep(1,ncol(x)),scale=FALSE,nbasis,basis.id,parallel=FALSE,cpus=1)
}


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension. 
         The range of input variables is scaled to [0,1].}
\item{time}{observed event time.}
\item{status}{status indicator. 1: dead; 0: alive.}
\item{wt}{weights for predictors. Default is \code{rep(1,ncol(x))}.}
\item{scale}{if \code{TRUE}, each predictor variable is rescaled to [0,1] interval. Dafault is \code{FALSE}.}
\item{nbasis}{number of "knots" to be used. Ignored when \code{basis.id} is specified.}
\item{basis.id}{index of observations to be used as "knots".}
\item{parallel}{parallelize task using \code{snowfall} package? Default is \code{FALSE}. Recommended when either sample size or number of predictor is large.}
\item{cpus}{number of available processor units. Default is \code{1}.}
}

\details{
In a proportional hazards model, the log-relative risk function is first assumed to have an additive form
\deqn{\lambda(t|x)=\lambda(t)\exp[\eta(x)]\\,~\eta(x)=\sum_{j=1}^p\eta_j(x_j),}


and then estimated via the penalized partial likelihood:
\deqn{-partial~likelihood/nobs+\lambda_0\sum_{j=1}^p\theta^{-1}_jw_j^2||\eta_j||^2, s.t.~\sum_{j=1}^p\theta_j<M.}

For large data sets, we can reduce the computational load of the optimizing problem by
selecting a subset of the observations of size \emph{nbais} as "knots", which reduces the dimension of the
kernel matrices from \emph{nobs} to \emph{nbasis}.
Unless specified via \code{basis.id} or \code{nbasis}, the number of "knots" is determined by \code{max}(35,12\emph{nobs}^{2/9}).


The weights can be specified based on either user's own discretion or adaptively computed from initial 
function estimates. See Storlie et al. (2011) for more discussions. One possible choice is to specify the weights
as the inverse \eqn{L_2} norm of initial function estimator, see \code{\link{SSANOVAwt.cox}}.
}

\value{
An object with S3 class "cosso".
\item{family}{type of regression model.}
\item{time}{observed event time.}
\item{status}{status indicator.}
\item{wt}{weights.}
\item{basis.id}{pre-specified or randomly selected index that are used as "knots".}
\item{tune}{a list containing tuning result}
}

\author{
Hao Helen Zhang and Chen-Yen Lin
}

\references{
Leng, C. and Zhang, H. H. (2006) "Model selection in nonparametric hazard regression", Nonparametric Statistics, \bold{18}, 417--429.
}

\seealso{\code{\link{tune.cosso.cox}}, \code{\link{plot.cosso}}, \code{\link{predict.cosso}}
}

\examples{
data(veteran)
t0=proc.time()
## Use half of the observations for demonstration
set.seed(27695)
train.id <- sort(sample(1:nrow(veteran),ceiling(nrow(veteran)/2)))
cossoCox <- cosso.cox(veteran[train.id,5:7],veteran[train.id,1],veteran[train.id,2],nbasis=20)
print((proc.time()-t0)[3])

\dontrun{
## Parallel Computing
## Not recommended in this example
data(veteran)
t0=proc.time()
## Use half of the observations for demonstration
set.seed(27695)
train.id <- sort(sample(1:nrow(veteran),ceiling(nrow(veteran)/2)))
cossoCox <- cosso.cox(veteran[train.id,5:7],veteran[train.id,1],veteran[train.id,2],parallel=TRUE,cpus=2)
print((proc.time()-t0)[3])

## Adaptive COSSO-Cox
data(veteran)
adawt <- SSANOVAwt.cox(veteran[,5:7],veteran[,1],veteran[,2])
ada.CoxObj <- cosso.cox(veteran[,5:&],veteran[,1],veteran[,2],adawt,FALSE,parallel=FALSE)
}
}
