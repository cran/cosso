\name{plot.cosso}
\alias{plot.cosso}
\title{
Plot method for COSSO object
}

\description{
Plot \eqn{L_2} norm solution path or main effects of selected functional components 
}

\usage{ 
\method{plot}{cosso}(x,M,plottype =c("Path","Functionals"),eps=1e-10,...) 
}


\arguments{
\item{x}{a cosso object}
\item{M}{a smoothing parameter value. M should be taken between 0 and p. Arguement required when \code{plottype="Functionals"}}
\item{plottype}{either \code{Path} (default) or \code{Functionals}. The Path plot shows the \eqn{L_2} norm path for each functional component as a function of smoothing parameter M. 
                The Functional plot shows the estimated functional components for each input variable at a particular smoothing parameter M. Abbreviations are allowed.}
\item{eps}{an effective zero, default is \code{1e-10}}
\item{...}{additional arguments for plot generic}
}


\value{
NULL
}

\references{
Lin, Y and Zhang, H. H. (2006). "Component Selection and Smoothing in Smoothing Spline Analysis of Variance Models," Annals of Statistics, \bold{34}, 2272--2297.

Storlie, C. B., Bondell, H. D., Reich, B. J. and Zhang, H. H. (2011). "Surface Estimation, Variable Selection, and the Nonparametric Oracle Property," Statistica Sinica, \bold{21}, 679--705.
}

\author{
Hao Helen Zhang  \email{hzhang@stat.ncsu.edu} }

\seealso{ \code{\link{predict.cosso}}
}



\examples{ 
data(ozone)
cossoObj <- cosso(x=ozone[,-1],y=ozone[,1],nbasis=100)
plot.cosso(cossoObj,plottype="Path")
plot.cosso(cossoObj,M=2,plottype="F")
}
