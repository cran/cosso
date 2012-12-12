\name{tune.cosso}
\alias{tune.cosso}
\title{
Compute K-fold cross-validated error or BIC score for COSSO
}

\description{
Compute the K-fold cross-validated mean squared prediction error or BIC score for COSSO.
}

\usage{ 
tune.cosso(object,type=c("BIC","CV"),folds=5,plot.it=TRUE) 
}


\arguments{
\item{object}{a cosso object}
\item{type}{type of tuning method}
\item{folds}{number of folds for corss-validation. Arguement required when \code{type="CV"}. Default is \code{5}}
\item{plot.it}{if \code{TRUE}, plot the cross-validated error or the BIC score curve}
}


\value{
\item{OptM}{the optimal smoothing parameter for M}
\item{OptLam}{the optimal smoothing parameter for \eqn{\lambda_0}}
\item{Mgrid}{a grid points for smoothing parameter M at which CV/BIC is computed}
\item{IC}{a list containing either cross-validated error or BIC score}
\item{L2norm}{functional \eqn{L_2}-norm computed at each grid point}
}

\references{
Lin, Y and Zhang, H. H. (2006). "Component Selection and Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.

Storlie, C. B., Bondell, H. D., Reich, B. J. and Zhang, H. H. (2011). "Surface estimation, variable selection, and the nonparametric oracle property", Statistica Sinica, \bold{21}, 679--705.
}

\author{
Hao Helen Zhang and Chen-Yen Lin}

\seealso{\code{\link{cosso}}, \code{\link{predict.cosso}}
}



\examples{ 
data(ozone)
cossoObj <- cosso(x=ozone[,2:5],y=ozone[,1],nbasis=30)
tune.cosso(cossoObj,type="BIC")

\dontrun{
## Tune by 5-fold CV
tune.cosso(cossoObj,type="CV")
}
}
