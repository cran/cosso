\name{summary.cosso}
\alias{summary.cosso}
\title{
Summary method for COSSO object}

\description{
Return a summary report for the fitted COSSO model
}

\usage{ 
\method{summary}{cosso}(object,eps=1e-10,...) 
}


\arguments{
\item{object}{object returned by cosso function}
\item{eps}{An effective zero, default is \code{1e-10}}
\item{...}{Additional arguments for summary generic}
}


\value{
\item{Rss}{mean of the squared residuals on the training data}
\item{Varlist}{list of selected function components}
\item{Df}{degree of freedom, trace of the smoothing matrix}
}

\references{ Lin, Y and Zhang, H.H. (2006) "Component Selection and
Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.}

\author{
Hao Helen Zhang \email{hzhang@stat.ncsu.edu} }

\seealso{ \code{\link{cosso}},  \code{\link{plot.cosso}},  \code{\link{predict.cosso}}}

\examples{
data(ozone) 
set.seed(2010)
# Randomly select one-third of the data as training data.
train_id <- sample(1:nrow(ozone),round(nrow(ozone)/3))
ozone_cosso <- cosso(x=ozone[train_id,-1],y=ozone[train_id,1],type="BIC")
summary.cosso(ozone_cosso)
}
