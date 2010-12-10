\name{predict.cosso} 
\alias{predict.cosso}

\title{
Make predictions or extract coefficients from a fitted COSSO model}

\description{
Make prediction for future observations
}

\usage{ 
\method{predict}{cosso}(object,xnew,type=c("fit","coefficients"),plot.it=TRUE,...) 
}


\arguments{
\item{object}{object returned by the cosso function}
\item{xnew}{future input values for response prediction}
\item{type}{If \code{type="fit"}, fitted values will be returned. If \code{type="coefficients"}, model coefficients will be returned. Abbreviations are allowed.}
\item{plot.it}{If \code{TRUE}, the fitted curve will be produced for each selected component in the final model.}
\item{...}{Additional arguments for predict generic}
}


\value{
Either a vector of fitted values, or a list containing estimated model coefficients. 
}

\references{ Lin, Y and Zhang, H.H. (2006) "Component Selection and
Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.}

\author{
Hao Helen Zhang \email{hzhang@stat.ncsu.edu} }

\seealso{ \code{\link{cosso}},  \code{\link{plot.cosso}}, \code{\link{summary.cosso}}
}

\examples{ 
data(ozone)
# Randomly select one-third of the data as training data.
set.seed(2010)
train_id <- sample(1:nrow(ozone),round(nrow(ozone)/3))
ozone_cosso <- cosso(x=ozone[train_id,-1],y=ozone[train_id,1],type="CV")
predict.cosso(ozone_cosso, ozone[-train_id,-1],type="fit")$value
}
