\name{plot.cosso}
\alias{plot.cosso}
\title{
Plot method for COSSO object
}

\description{
Plot the fitted curve for each selected function component in the final model.
}

\usage{ 
\method{plot}{cosso}(x,eps=1e-10,...) 
}


\arguments{
\item{x}{object returned by cosso function}
\item{eps}{An effective zero, default is \code{1e-10}}
\item{...}{Additional arguments for plot generic}
}


\value{
NULL
}

\references{ Lin, Y and Zhang, H.H. (2006) "Component Selection and
Smoothing in Smoothing Spline Analysis of Variance Models", Annals of Statistics, \bold{34}, 2272--2297.}

\author{
Hao Helen Zhang  \email{hzhang@stat.ncsu.edu} }

\seealso{\code{\link{cosso}}, \code{\link{predict.cosso}}, \code{\link{summary.cosso}}
}


\examples{ 
data(ozone)
set.seed(2010)
# Randomly select one-third of the data as training data.
train_id <- sample(1:nrow(ozone),round(nrow(ozone)/3))
ozone_cosso <- cosso(x=ozone[train_id,-1],y=ozone[train_id,1],type="BIC")
plot.cosso(ozone_cosso)
}
