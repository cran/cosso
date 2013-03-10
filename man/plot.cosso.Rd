\name{plot.cosso}
\alias{plot.cosso}
\title{
Plot method for COSSO object
}

\description{
Plot \eqn{L_2} norm solution path or main effects of selected functional components 
}

\usage{ 
\method{plot}{cosso}(x,M,plottype =c("Path","Functionals"),eps=1e-7,...) 
}


\arguments{
\item{x}{a cosso object.}
\item{M}{a smoothing parameter value. Arguement required when \code{plottype="Functionals"}.}
\item{plottype}{either \code{Path} (default) or \code{Functionals}. The Path plot shows the \eqn{L_2} norm path for each functional component as a function of smoothing parameter M. 
                The Functional plot shows the estimated functional components for each input variable at a particular smoothing parameter M. Abbreviations are allowed.}
\item{eps}{an effective zero, default is \code{1e-7}.}
\item{...}{additional arguments for plot generic.}
}


\value{
NULL
}

\author{
Hao Helen Zhang and Chen-Yen Lin}

\seealso{ \code{\link{predict.cosso}}
}



\examples{ 
set.seed(20130310)
x=cbind(rbinom(200,1,.7),matrix(runif(200*7,0,1),nc=7))
trueProb=1/(1+exp(-x[,1]-sin(2*pi*x[,2])-5*(x[,4]-0.4)^2))
y=rbinom(200,1,trueProb)

B.Obj=cosso(x,y,family="Bin") 
plot.cosso(B.Obj,plottype="Path")
plot.cosso(B.Obj,M=2,plottype="Func")
}
