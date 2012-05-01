\name{cosso.qr} 
\alias{cosso.qr}
\title{
Fit a nonparaletric quantile regression model with COSSO regularization
}

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
The weights can be specified based on either user's own discretion or adaptively computed from initial 
function estimates. See Storlie et al. (2011) for more discussions. One possible choice is to specify the weights
as the inverse \eqn{L_2} norm of initial function estimator, see \code{\link{KQRwt}}.
}

\value{
An object with S3 class "cosso".
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
set.seed(27695)
## Use one third observations as training set
train_id   <- sample(1:nrow(ozone),round(nrow(ozone)/3))
cossoqrObj <- cosso.qr(x=ozone[train_id,-1],y=ozone[train_id,1],tau=0.5)
plot.cosso(cossoqrObj,plottype="Path")

## Adaptive COSSO-QR
adaptive.wt <- KQRwt(ozone[train_id,-1],ozone[train_id,1],0.5)
acossoqrObj <- cosso.qr(x=ozone[train_id,-1],y=ozone[train_id,1],tau=0.5,wt=adaptive.wt)
}
