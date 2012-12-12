\name{tune.cosso.qr}
\alias{tune.cosso.qr}
\title{
Compute K-fold cross-validated error for COSSO quantile regression
}

\description{
Compute K-fold cross-validated check error for COSSO quantile regression model.
}

\usage{
tune.cosso.qr(object,folds=5,plot.it=TRUE,parallel=FALSE,cpus=1)
}


\arguments{
\item{object}{a cosso object.}
\item{folds}{number of folds. Default is \code{5}.}
\item{plot.it}{if \code{TRUE}, plot the cross-validation error curve.}
\item{parallel}{parallelize task using \code{snowfall} package? Default is \code{FALSE}. Recommended when either sample size or number of predictors is large.}
\item{cpus}{number of available cpu unit. Default is \code{1}. Arguement required when parallel=\code{TRUE}.}
}


\value{
\item{OptM}{the selected smoothing parameter for M.}
\item{OptLam}{the selected smoothing parameter for \eqn{\lambda_0}.}
\item{Mgrid}{a grid points for smoothing parameter M at which cross-validation error is computed.}
\item{IC}{a vector containing cross-validation error.}
\item{L2norm}{functional \eqn{L_2}-norm for each input variable.}
}


\author{
Hao Helen Zhang and Chen-Yen Lin}

\seealso{\code{\link{cosso.qr}}, \code{\link{predict.cosso}}
}



\examples{
data(ozone)
t0=proc.time()
set.seed(27695)
## Use one third observations for demonstration
train_id <- sample(1:nrow(ozone),round(nrow(ozone)/3))
cossoqrObj <- cosso.qr(x=ozone[train_id,2:5],y=ozone[train_id,1],tau=0.5)
tune.cosso.qr(cossoqrObj)
(proc.time()-t0)[3]

\dontrun{
## Parallel Computing
tune.cosso.qr(cossoqrObj,parallel=TRUE,cpus=2)
}
}
