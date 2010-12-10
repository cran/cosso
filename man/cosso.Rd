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
cosso(x,y,wt=rep(1,ncol(x)),type=c("BIC","CV"),folds=5)
}


\arguments{
\item{x}{input matrix; the number of rows is sample size, the number of columns is the data dimension. 
         The range of input variables is scaled to [0,1].}
\item{y}{response vector}
\item{wt}{weights for predictors. Default is \code{rep(1,ncol(x))}}
\item{type}{tuning criterion}
\item{folds}{the number of folds for cross-validation. Argument required when \code{type="CV"}. Default is \code{5}}
}

\details{
The weights can be specified based on either user's own discretion or adaptively computed from initial 
function estimates. See Storlie et al. (2011) for more discussions. One possible choice is to specify the weights
as the inverse \eqn{L_2} norm of initial function estimator, see \code{\link{compwt}}.

}

\value{
\item{b}{the estimated intercept}
\item{cc}{the estimated coefficients for kernel representers}
\item{Df}{degree of freedom, the trace of the smoothing matrix}
\item{theta}{the estimated scale parameter for each component}
\item{wt}{pre-specified weights}
\item{x}{the input matrix}
\item{y}{the response vector}
}

\references{
Lin, Y and Zhang, H. H. (2006) "Component Selection and Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.
}

\author{
Hao Helen Zhang \email{hzhang@stat.ncsu.edu}
}

\seealso{ \code{\link{plot.cosso}}, \code{\link{predict.cosso}}, \code{\link{summary.cosso}}
}

\examples{
data(ozone)
\dontrun{## Fit cosso with 5-CV tuning
cossoCV <- cosso(x=ozone[,-1],y=ozone[,1],type="CV",folds=5)
summary.cosso(cossoCV)
plot.cosso(cossoCV)}

## Fit cosso with BIC tuning
set.seed(2010)
# Randomly select one-third of the data as training data.
train_id <- sample(1:nrow(ozone),round(nrow(ozone)/3))
cossoBIC <- cosso(x=ozone[train_id,-1],y=ozone[train_id,1],type="BIC")
summary.cosso(cossoBIC)
plot.cosso(cossoBIC)

## Fit adaptive cosso
wt <- compwt(ozone[train_id,-1],ozone[train_id,1],rep(1,ncol(ozone)-1),1)
acosso <- cosso(x=ozone[train_id,-1],y=ozone[train_id,1],wt=wt,type="BIC")
summary.cosso(acosso)
}
