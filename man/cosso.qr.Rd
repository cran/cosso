\name{cosso.qr} 
\alias{cosso.qr}
\title{
Estimate the quantile function by Smoothing Splines with COSSO penalty}

\description{
Fit a penalized quantile regression model with COSSO regularization. The regularization path is computed
for the COSSO penalty at a grid of values for the regularization parameter M.}

\usage{
cosso.qr(x,y,tau,wt=rep(1,ncol(x)),scale=FALSE,parallel=FALSE,cpus=1)
}


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension. 
         The range of input variables is scaled to [0,1].}
\item{y}{response vector.}
\item{tau}{the quantile to be estimated, a number strictly between 0 and 1.}
\item{wt}{weights for predictors. Default is \code{rep(1,ncol(x))}.}
\item{scale}{if \code{TRUE}, each predictor variable is rescaled to [0,1] interval. Dafault is \code{FALSE}.}
\item{parallel}{parallelize task using \code{snowfall} package? Default is \code{FALSE}. Recommended when sample size is large.}
\item{cpus}{number of available processor units. Default is \code{1}.}
}

\details{
The quantile function is first assumed to have an additive form:
\deqn{\eta(x)=\sum_{j=1}^p\eta_j(x_j),}

then estimated by minimizing the objective function:
\deqn{check~loss/nobs+\lambda_0\sum_{j=1}^p\theta^{-1}_jw_j^2||\eta_j||^2, s.t. \sum_{j=1}^p\theta_j<M.}


The weights can be specified based on either user's own discretion or adaptively computed from initial 
function estimates. See Storlie et al. (2011) for more discussions. One possible choice is to specify the weights
as the inverse \eqn{L_2} norm of initial function estimator, see \code{\link{KQRwt}}.
}

\value{
An object with S3 class "cosso".
\item{family}{type of regression model.}
\item{tau}{quantile to be estimated.}
\item{wt}{weights.}
\item{x}{input matrix.}
\item{y}{response vector.}
\item{basis.id}{Indices of observations used to compute the kernel matrices. All observations will be used in quantile regression.}
\item{tune}{a list containing tuning result}
}

\author{
Hao Helen Zhang and Chen-Yen Lin
}

\seealso{\code{\link{tune.cosso.qr}}, \code{\link{plot.cosso}}, \code{\link{predict.cosso}}
}

\examples{
data(ozone)
t0 <- proc.time()
set.seed(27695)
## Use one third observations for demonstration
train_id   <- sample(1:nrow(ozone),round(nrow(ozone)/3))
cossoqrObj <- cosso.qr(x=ozone[train_id,2:5],y=ozone[train_id,1],tau=0.5)
print((proc.time()-t0)[3])

\dontrun{
data(ozone)
## Parallel Computing
t0 <- proc.time()
set.seed(27695)
train_id   <- sample(1:nrow(ozone),round(nrow(ozone)/3))
cossoqrObj <- cosso.qr(x=ozone[train_id,2:9],y=ozone[train_id,1],tau=0.5,parallel=TRUE,cpus=2)
print((proc.time()-t0)[3])

data(ozone)
set.seed(27695)
train_id   <- sample(1:nrow(ozone),round(nrow(ozone)/3))
## Adaptive COSSO-QR
adaptive.wt <- KQRwt(ozone[train_id,2:9],ozone[train_id,1],0.5)
acossoqrObj <- cosso.qr(x=ozone[train_id,2:9],y=ozone[train_id,1],tau=0.5,wt=adaptive.wt)
}
}
