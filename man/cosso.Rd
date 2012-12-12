\name{cosso} 
\alias{cosso}
\title{
Estimate the mean regression function for Gaussian response using Smmohting Spines with COSSO penalty
}

\description{
Fit COSSO and adaptive COSSO models for Gaussian response.
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
The mean regression function is first assumed to have an additive form
\deqn{\eta(x)=\sum_{j=1}^p\eta_j(x_j),}

then estimated by minimizing the objective function:
\deqn{RSS/nobs+\lambda_0\sum_{j=1}^p\theta^{-1}_jw_j^2||\eta_j||^2, s.t.~\sum_{j=1}^p\theta_j<M.}


For large data sets, we can reduce the computational load of the optimizing problem by
selecting a subset of the observations of size \emph{nbais} as "knots", which reduces the dimension of the
kernel matrices from \emph{nobs} to \emph{nbasis}.
Unless specified via \code{basis.id} or \code{nbasis}, the default number of "knots" is the sample size (\emph{nobs}).


The weights can be specified based on either user's own discretion or adaptively computed from initial 
function estimates. See Storlie et al. (2011) for more discussions. One possible choice is to specify the weights
as the inverse \eqn{L_2} norm of initial function estimator, see \code{\link{SSANOVAwt}}.
}

\value{
An object with S3 class "cosso".
\item{family}{type of regression model.}
\item{x}{the input matrix}
\item{y}{the response vector}
\item{Kmat}{an array containing kernel matrices for each input variables}
\item{basis.id}{Indices of observations used as "knots"}
\item{wt}{weights}
\item{tune}{a list containing prelminary tuning result}
}

\references{
Lin, Y. and Zhang, H. H. (2006) "Component Selection and Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.

Storlie, C. B., Bondell, H. D., Reich, B. J. and Zhang, H. H. (2011) "Surface Estimation, Variable Selection, and the Nonparametric Oracle Property", Statistica Sinica, \bold{21}, 679--705.
}

\author{
Hao Helen Zhang}

\seealso{ \code{\link{plot.cosso}}, \code{\link{predict.cosso}}, \code{\link{tune.cosso}}
}

\examples{
data(ozone)
## Fit cosso
## Use 50 observations as knots
t0=proc.time()
## Use half of the observations for demonstration
set.seed(27695)
train.id <- sort(sample(1:nrow(ozone),ceiling(nrow(ozone)/2)))
cossoObj <- cosso(x=ozone[train.id,2:5],y=ozone[train.id,1],nbasis=50)
print((proc.time()-t0)[3])

\dontrun{
## Use all observations as knots
t0=proc.time()
## Use half of the observations for demonstration
set.seed(27695)
train.id <- sort(sample(1:nrow(ozone),ceiling(nrow(ozone)/2)))
cossoObj <- cosso(x=ozone[train.id,2:5],y=ozone[train.id,1])
print((proc.time()-t0)[3])


## Fit adaptive cosso
adaptive.wt <- SSANOVAwt(ozone[,-1],ozone[,1])
acossoObj <- cosso(x=ozone[,-1],y=ozone[,1],wt=adaptive.wt,nbasis=ceiling(nrow(ozone)/5))
}
}
