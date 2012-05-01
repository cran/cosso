\name{cosso} 
\alias{cosso}
\title{
Fit the additive COSSO and adaptive COSSO models.
}

\description{
COSSO is a regularization method for variable selection and function estimation in multivariate
nonparametric regression models. By imposing a soft-thresholding type penalty onto function components, 
the COSSO solution is sparse and hence able to identify important variables. 
The method is developed in the framework of smoothing spline ANOVA.}

\usage{
cosso(x,y,wt=rep(1,ncol(x)),scale=FALSE,nbasis,basis.id,n.step=2*ncol(x))
}


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension. 
         The range of input variables is scaled to [0,1].}
\item{y}{response vector}
\item{wt}{weights for predictors. Default is \code{rep(1,ncol(x))}}
\item{scale}{if TRUE, each predictor variable is rescaled to [0,1] interval. Dafault is \code{FALSE}.}
\item{basis.id}{index designating selected "knots".}
\item{nbasis}{number of "knots" to be selected. Ignored when \code{basis.id} is provided.}
\item{n.step}{maximum iteration number in fiding solution path.}
}

\details{
The weights can be specified based on either user's own discretion or adaptively computed from initial 
function estimates. See Storlie et al. (2011) for more discussions. One possible choice is to specify the weights
as the inverse \eqn{L_2} norm of initial function estimator, see \code{\link{SSANOVAwt}}.

A subset of the observations can be selected as "knots" to faciliate computation. Default number of "knots" is the sample size.

}

\value{
\item{x}{the input matrix}
\item{y}{the response vector}
\item{Kmat}{an array containing kernel matrices for each input variables}
\item{basis.id}{Indices of observations used as "knots"}
\item{wt}{pre-specified weights}
\item{tune}{a list containing prelminary tuning result}
}

\references{
Lin, Y. and Zhang, H. H. (2006) "Component Selection and Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.

Zhang, H. H. and Lin, Y. (2006) "Component Selection and Smoothing for Nonparametric Regression in Exponential Families," Statistica Sinica, \bold{16}, 1021--1041.

Storlie, C. B., Bondell, H. D., Reich, B. J. and Zhang, H. H. (2011) "Surface Estimation, Variable Selection, and the Nonparametric Oracle Property", Statistica Sinica, \bold{21}, 679--705.
}

\author{
Hao Helen Zhang \email{hzhang@stat.ncsu.edu}
}

\seealso{ \code{\link{plot.cosso}}, \code{\link{predict.cosso}}, \code{\link{tune.cosso}}
}

\examples{
data(ozone)
## Fit cosso
## Use one third observations as knots
cossoObj <- cosso(x=ozone[,-1],y=ozone[,1],nbasis=ceiling(nrow(ozone)/3))
plot.cosso(cossoObj,plottype="Path")

## Fit adaptive cosso
adaptive.wt <- SSANOVAwt(ozone[,-1],ozone[,1])
acossoObj <- cosso(x=ozone[,-1],y=ozone[,1],wt=adaptive.wt,nbasis=ceiling(nrow(ozone)/3))
}
