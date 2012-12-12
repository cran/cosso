\name{predict.cosso} 
\alias{predict.cosso}

\title{
Make predictions or extract coefficients from a COSSO model}

\description{
Make prediction for future observations or extract the model parameters at a particular smoothing parameter.
}

\usage{ 
\method{predict}{cosso}(object,xnew,M,type=c("fit","coefficients"),...) 
}


\arguments{
\item{object}{a cosso object.}
\item{xnew}{future input values for response prediction}
\item{M}{a smoothing parameter value. M should be taken between 0 and p. If not provided, a cross-validation procedure will be carried out to select an appropriate value.}
\item{type}{if \code{type="fit"} (default), fitted values will be returned. If \code{type="coefficients"}, model coefficients will be returned. Abbreviations are allowed.}
\item{...}{additional arguments for predict function.}
}


\value{
Either a vector of fitted values, or a list containing estimated model coefficients. 

When \code{type="coefficients"}, three sets of coefficients will be returned.
\item{Intercept}{the estimated intercept}
\item{coefs}{the estimated coefficients for kernel representers}
\item{theta}{the estimated scale parameters for each functional component}
}


\author{
Hao Helen Zhang and Chen-Yen Lin}

\seealso{\code{\link{plot.cosso}}
}

\examples{ 
data(ozone)
## Use half of the observations for demonstration
set.seed(27695)
train.id <- sort(sample(1:nrow(ozone),ceiling(nrow(ozone)/2)))
cossoObj <- cosso(x=ozone[train.id,2:5],y=ozone[train.id,1],nbasis=50)
## Compute model parameters
predict.cosso(cossoObj,M=2,type="c")
## Make prediction on full data
predict.cosso(cossoObj,ozone[,2:5],M=2,type="c")
}
