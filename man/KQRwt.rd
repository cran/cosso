\name{KQRwt}
\alias{KQRwt}
\title{
Compute adaptive weights for Quantile regression
}

\description{
A preliminary Kernel Quantile Regression model \eqn{\tilde{\eta}} is first obtained by
minimizing a Kernel Quantile Regression type of objective function, 
and then the weight for the j-th function component is computed as \eqn{||P_j \tilde{\eta}||^{-1}}. 
Here \eqn{P_j} denotes the projection operator to the subspace.
}

\usage{ KQRwt(x,y,tau,cand.lam0,folds=5,parallel=FALSE,cpus=1) 
}


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension. 
         The range of input variables is scaled to [0,1].}
\item{y}{response vector}
\item{tau}{the quantile to be estimated, a number strictly between 0 and 1.}
\item{cand.lam0}{a list of tuning parameter \eqn{\lambda0}. Default is 2^seq(-14,-20,-0.75)}
\item{folds}{number of folds for corss-validation. Default is \code{5}.}
\item{parallel}{parallelize task using \code{snowfall} package? Default is \code{FALSE}. Recommended when either sample size or number of predictors is large.}
\item{cpus}{number of available cpu unit. Default is \code{1}. Arguement required when parallel=\code{TRUE}.}
}

\details{
The initial quantile function is estimated by the objective function:
\deqn{check~loss/nobs+\lambda_0\sum_{j=1}^p||P_j\eta||^2}.
}



\value{
\item{wt}{The adaptive weights}
}


\author{
Hao Helen Zhang and Chen-Yen Lin}


\references{
Y. Li, Liu, Y. and Zhu, J. (2007) "Quantile regression in reproducing kernel Hilbert spaces," Journal of the American Statistical Association, \bold{477}, 255--267.
}


\seealso{ \code{\link{cosso.qr}} }

\examples{
## Parallel Computing
data(ozone) 
set.seed(27695)
## Use half observations as training set
train_id <- sample(1:nrow(ozone),round(nrow(ozone)/2))
QRwt     <- KQRwt(x=ozone[train_id,-1],y=ozone[train_id,1],tau=0.5)
print(QRwt)

\dontrun{
## Use all observations and conduct parallel computing
QRwt <- KQRwt(ozone[,-1],ozone[,1],0.5,parallel=TRUE,cpus=3)
print(QRwt)
}
}
