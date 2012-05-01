\name{KQRwt}
\alias{KQRwt}
\title{
Compute initial weights used in the penalty term of the adaptive COSSO quantile regression model}

\description{
A preliminary Kernel Quantile Regression model \eqn{f} is first obtained, 
and then the weight for the jth function component is computed as \eqn{||P_j f||^{-\gamma}}. 
Here \eqn{P_j} denotes the projection operator to the subspace.
}

\usage{ KQRwt(x,y,tau,gampow=1,folds=5,parallel=FALSE,cpus=1) }


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension. 
         The range of input variables is scaled to [0,1].}
\item{y}{response vector}
\item{tau}{the quantile to be estimated, a number strictly between 0 and 1}
\item{gampow}{power of the \eqn{L_2} norm. Default is \code{1}}
\item{folds}{number of folds for corss-validation. Default is \code{5}}
\item{parallel}{parallelize task using \code{snowfall} package? Default is \code{FALSE}. Recommended when sample size is large.}
\item{cpus}{number of available cpu unit. Default is \code{1}}
}


\value{
\item{wt}{Reciprocal of \eqn{L_2} norm. The adaptive weights can be used in the COSSO Quantile Regression}
}


\author{
Hao Helen Zhang }


\references{
Y. Li, Liu, Y. and Zhu, J. (2007) "Quantile regression in reproducing kernel Hilbert spaces," Journal of the American Statistical Association, \bold{477}, 255--267.
}


\seealso{ \code{\link{KQR}},\code{\link{cosso.qr}} }

\examples{
data(ozone) 
set.seed(27695)
t0 <- proc.time()
QRwt <- KQRwt(ozone[,-1],ozone[,1],0.5)
(proc.time()-t0)[3]
print(QRwt)

## Parallel Computing
set.seed(27695)
t0<-proc.time()
QRwt <- KQRwt(ozone[,-1],ozone[,1],0.5,parallel=TRUE,cpus=2)
(proc.time()-t0)[3]
print(QRwt)
}
