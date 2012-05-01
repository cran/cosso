\name{KQR}
\alias{KQR}
\title{Kernel Quantile Regression Method}

\description{
Build a Kernel Quantile Regression Model 
}

\usage{ KQR(x,y,tau,lam0s,folds=5,parallel=FALSE,cpus=1) }


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension. 
         The range of input variables is scaled to [0,1].}
\item{y}{response vector.}
\item{tau}{the quantile to be estimated, a number strictly between 0 and 1.}
\item{lam0s}{tuning grid points for smoothing parameter \eqn{\lambda}. The default value is \code{2^seq(-14,-20,by=-0.5)}.}
\item{folds}{number of folds for corss-validation. Default is \code{5}.}
\item{parallel}{parallelize task using \code{snowfall} package? Default is \code{FALSE}. Recommended when sample size is large.}
\item{cpus}{number of available cpu unit. Default is \code{1}. Arguement required when parallel=\code{TRUE}.}
}


\value{
\item{coefs}{estimated coefficients for kernel representers}
\item{intercept}{estimated intercept}
\item{quantile}{estimated quantile function value at design points.}
\item{L2norm}{\eqn{L_2} norm for each of the functional component.}
\item{OptLam}{selected smoothing parameter.}
}


\author{
Hao Helen Zhang and Chen-Yen Lin
}


\references{
Y. Li, Liu, Y. and Zhu, J. (2007) "Quantile regression in reproducing kernel Hilbert spaces," Journal of the American Statistical Association, \bold{477}, 255--267.
}


\seealso{ \code{\link{KQRwt}} }

\examples{
data(ozone) 
set.seed(27695)
kqrObj <- KQR(ozone[,-1],ozone[,1],0.4)
kqrObj$quantile

## Parallel Computing
set.seed(27695)
kqrObj.p <- KQR(ozone[,-1],ozone[,1],0.4,parallel=TRUE,cpus=2)
kqrObj.p$quantile
}
