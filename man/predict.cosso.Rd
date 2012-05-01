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
\item{object}{a cosso object}
\item{xnew}{future input values for response prediction}
\item{M}{a smoothing parameter value. M should be taken between 0 and p. If not provided, a cross-validation procedure will be carried out to select an appropriate value.}
\item{type}{if \code{type="fit"} (default), fitted values will be returned. If \code{type="coefficients"}, model coefficients will be returned. Abbreviations are allowed.}
\item{...}{additional arguments for predict generic}
}


\value{
Either a vector of fitted values, or a list containing estimated model coefficients. 

When \code{type="coefficients"}, three sets of coefficients will be returned.
\item{Intercept}{the estimated intercept}
\item{coefs}{the estimated coefficients for kernel representers}
\item{theta}{the estimated scale parameters for each functional component}
}

\references{
Lin, Y and Zhang, H. H. (2006). "Component Selection and Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.

Storlie, C. B., Bondell, H. D., Reich, B. J. and Zhang, H. H. (2011). "Surface estimation, variable selection, and the nonparametric oracle property", Statistica Sinica, \bold{21}, 679--705.
}

\author{
Hao Helen Zhang \email{hzhang@stat.ncsu.edu} }

\seealso{\code{\link{plot.cosso}}
}

\examples{ 
data(ozone)
ozone_cosso <- cosso(x=ozone[,-1],y=ozone[,1],nbasis=100)
predict.cosso(ozone_cosso,ozone[,-1],M=2,type="fit")
predict.cosso(ozone_cosso,M=2,type="c")
}
