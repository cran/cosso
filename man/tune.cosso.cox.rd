\name{tune.cosso.cox}
\alias{tune.cosso.cox}
\title{
Compute the approximate cross-validated score for COSSO-Cox model.
}

\description{
Compute the approximate cross-validated score (Leng and Zhang, 2006) for COSSO-Cox model.
}

\usage{
tune.cosso.cox(object,plot.it=TRUE,parallel=FALSE,cpus=1)
}


\arguments{
\item{object}{a cosso object.}
\item{plot.it}{if \code{TRUE}, plot the approximate cross-validation score against smoothing parameter M.}
\item{parallel}{parallelize task using \code{snowfall} package? Default is \code{FALSE}. Recommended when either sample size or number of predictors is large.}
\item{cpus}{number of available cpu unit. Default is \code{1}. Arguement required when parallel=\code{TRUE}.}
}


\value{
\item{OptM}{the selected smoothing parameter for M.}
\item{OptLam}{the selected smoothing parameter for \eqn{\lambda_0}.}
\item{Mgrid}{a grid points for smoothing parameter M at which cross-validation error is computed.}
\item{ACV}{a vector containing approximate cross-validation score.}
\item{L2norm}{functional \eqn{L_2}-norm for each input variable.}
}


\author{
Hao Helen Zhang and Chen-Yen Lin}


\references{
Leng, C. and Zhang, H. H. (2006) "Model selection in nonparametric hazard regression", Nonparametric Statistics, \bold{18}, 417--429.
}


\seealso{\code{\link{cosso.cox}}, \code{\link{predict.cosso}}
}



\examples{
data(veteran)
t0 <- proc.time()
## Use half of the observations for demonstration
set.seed(27695)
train.id <- sort(sample(1:nrow(veteran),ceiling(nrow(veteran)/2)))
cossoCox <- cosso.cox(x=veteran[train.id,5:7],time=veteran[train.id,1],status=veteran[train.id,2],nbasis=30)
tuneOnj <- tune.cosso.cox(cossoCox,plot.it=TRUE)
print((proc.time()-t0)[3])

\dontrun{
## Parallel Computing
## Not recommended in this example
t0 <- proc.time()
cossoCox <- cosso.cox(x=veteran[,5:7],time=veteran[,1],status=veteran[,2],nbasis=25)
tune.cosso.cox(cossoqrObj,parallel=TRUE,cpus=2)
print((proc.time()-t0)[3])
}

}
