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
\item{x}{a cosso object}
\item{M}{a smoothing parameter value. M should be taken between 0 and p. Arguement required when \code{plottype="Functionals"}}
\item{plottype}{either \code{Path} (default) or \code{Functionals}. The Path plot shows the \eqn{L_2} norm path for each functional component as a function of smoothing parameter M. 
                The Functional plot shows the estimated functional components for each input variable at a particular smoothing parameter M. Abbreviations are allowed.}
\item{eps}{an effective zero, default is \code{1e-7}}
\item{...}{additional arguments for plot generic}
}


\value{
NULL
}

\author{
Hao Helen Zhang and Chen-Yen Lin}

\seealso{ \code{\link{predict.cosso}}
}



\examples{ 
data(ozone)
## Use half of the observations for demonstration
set.seed(27695)
train.id <- sort(sample(1:nrow(ozone),ceiling(nrow(ozone)/2)))
cossoObj <- cosso(x=ozone[train.id,2:5],y=ozone[train.id,1],nbasis=50)
plot.cosso(cossoObj,plottype="Path")
plot.cosso(cossoObj,M=2,plottype="Functionals")
}
