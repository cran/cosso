.packageName <- "cosso"

##generate spline kernel (including linear terms)
genK <- function(x1,x2)
{
m = length(x1)
n = length(x2)
b = abs(x1%*%matrix(1,ncol=n)-matrix(1,nrow=m)%*%t(x2))
k1s = x1 - 0.5
k1t = x2 - 0.5
k2s = (k1s^2 - 1/12)/2
k2t = (k1t^2 - 1/12)/2
K = k1s%*%t(k1t)+k2s%*%t(k2t)-((b-0.5)^4-(b-0.5)^2/2+7/240)/24
}


gram <- function(x1,x2,mscale)
 {   n1 <- dim(as.matrix(x1))[1]
     d  <- length(mscale)
     n2 <- dim(as.matrix(x2))[1]
     KK <- matrix(0,n1,n2)
     if (d==1)
       {KK = mscale*genK(x1,x2)}
     else
       {for (j in 1:d)
          { KK = KK + mscale[j]*genK(x1[,j],x2[,j])}
       }
     return(KK)
 }

bigGram <- function(x1,x2)
{ n1 <- nrow(x1)
  n2 <- nrow(x2)
  d <- ncol(x1)
  Gram <- array(0,c(n1,n2,d))
  for(j in 1:d){
    Gram[,,j] <-  genK(x1[,j], x2[,j])
  }
  return(Gram)
}  

wsGram <- function(Gramat,mscale)
 {  n1 <- dim(Gramat)[1]
    n2 <- dim(Gramat)[2]
    d <- dim(Gramat)[3]
    KK <- matrix(0,n1,n2)
    for (j in 1:d)
    { KK = KK + mscale[j]*Gramat[,,j]}
    return(KK)
 }
 

##Main COSSO function
cosso <- function(x,y,wt=rep(1,ncol(x)),type=c("BIC","CV"), folds=5)
 {   n <- dim(as.matrix(x))[1]
     d <- dim(as.matrix(x))[2]
     type <- match.arg(type) 
     Gramat <- bigGram(x,x)
     if (type=="BIC") 
       {est=bicadd(Gramat,y,wt)}       
     if (type=="CV")
       {est=cvadd(Gramat,y,wt,folds)}
     z <- list(x=x,cc=rep(0,n),b=0,theta=rep(1,d),wt=rep(1,d),Df=d)
     z$x = x
     z$y = y
     z$cc = est$solution[1:n]
     z$b = est$solution[n+1]
     z$theta = est$solution[(n+2):(n+1+d)]
     z$wt = wt
     z$Df = est$Df
     class(z)="cosso"
     return(z)
 }

predict.cosso <- function(object, xnew, type=c("fit","coefficients"),plot.it=TRUE,...)
  { x <- object$x
    cc <- object$cc
    b <- object$b
    theta <- object$theta
    wt <- object$wt
    type=match.arg(type)
    if (type=="fit")
     {Knew = gram(xnew, x, theta/(wt^2))
      value = Knew%*%cc+b
      predictor = list(value=value)
     }
    if (type=="coefficients")
        predictor = list(coef=cc,b=b,theta=theta)
    if (plot.it)
        plot.cosso(object)
    class(predictor)="predict.cosso"       
    return(predictor)
  } 


summary.cosso <- function(object,eps=1e-10,...)
 {  x <- object$x
    y <- object$y
    cc <- object$cc
    b <- object$b
    theta <- object$theta
    d <- length(theta)
    wt <- object$wt
    Df <- object$Df
    n = dim(as.matrix(x))[1]
    fpred <- predict.cosso(object,x,type="fit",plot.it=FALSE)
    Rss <- mean((y-fpred$value)^2) 
    Varlist <- (1:d)[abs(theta)>eps]
    if (is.null(names(x))==0)
      Varlist <- names(x)[Varlist]
    result <- list(Rss=Rss,Df=Df,Varlist=Varlist)
    class(result)="summary.cosso"
    return(result)
}

plot.cosso <- function(x,eps=1e-10,...)
  { 
    object<-x
    x <- object$x
    cc <- object$cc
    b <- object$b
    wt <- object$wt
    theta <- object$theta
    d <- length(theta)
    sid <- (1:d)[abs(theta)>eps]
    fsize <- length(sid)
    n = dim(as.matrix(x))[1]
    fits <- matrix(0,n,d)
    for (jj in 1:d)
     {fits[,jj]=genK(x[,jj],x[,jj])%*%cc*theta[jj]/(wt[jj]^2)}
    frange = range(fits)
    plotrow <- ceiling(fsize/2)
    plotcol <- 2
    par(mfrow=c(plotrow,plotcol))
    for (jj in 1:fsize)  #only plot the selected components
      {sjj = sid[jj]
       jorder = order(x[,sjj])
       if (is.null(names(x))==0)
	plot(x[jorder,sjj],fits[jorder,sjj],type="l",xlim=c(0,1),ylim=frange,xlab=names(x)[sjj],ylab=paste("f(",names(x)[sjj],")",sep=""))
       else
       	plot(x[jorder,sjj],fits[jorder,sjj],type="l",xlim=c(0,1),ylim=frange,xlab=paste("x",eval(sjj),sep=""),ylab=paste("f(x",eval(sjj),")",sep=""))
       }
  }

## solve the standard smoothing spline
sspline <- function(Gramat,y,mscale,lam)
  { n <- length(y)
    Kmat = wsGram(Gramat,mscale)
    bigK = matrix(1,n+1,n+1)
    bigK[1:n,1:n] = Kmat+lam*diag(1,n)
    bigK[n+1,n+1]= 0
    cb = solve(bigK,c(y,0))
    Df = compdf(Kmat,y,lam)
    result=list(cb=cb,Df=Df) 
    return(result)
 }


## split n samples into K folds, 
cvsplit <- function(n,folds)
  { fsize <- floor(n/folds) #average size of each fold, the last fold is big
    splits = fsize*rep(1,folds)
    nextra = n-folds*fsize
    if (nextra>0) 
      {splits[1:nextra] = splits[1:nextra]+1}
    return(splits)
  }


## input: n, K
## output: the fold size 
cvlam <- function(Gramat,y,mscale,folds)
  { n <- length(y)
    d <- dim(Gramat)[3]
    splits = cvsplit(n,folds) 
    pp = sample(n,replace=FALSE) 
    LL <- 10
    lam0 <- 2^(-(1:LL))
    cverr <- matrix(0,LL,folds)
    for (kk in 1:folds)
    { tmp = 0*(kk==1)+sum(splits[1:(kk-1)])*(kk>1) #starting index 
      ntrain = n-splits[kk] #splits[kk] is the tuning size, ntrain = training size
      tuneid = pp[(tmp+1):(tmp+splits[kk])]
      ytune = y[tuneid]
      ytrain = y[-tuneid]
      trainGramat = Gramat[-tuneid,-tuneid,1:d]
      predGramat = Gramat[tuneid,-tuneid,1:d]
      Kpred = wsGram(predGramat,mscale)
      for (jj in 1:LL)
        {cb = sspline(trainGramat,ytrain,mscale,lam0[jj])$cb 
         predict = Kpred%*%cb[1:ntrain]+cb[ntrain+1]
         cverr[jj,kk] = sum((ytune-predict)^2)
        }
    }
      meancv = apply(cverr,1,sum)/n #the total leave-out error/n
      return(lam0[which.min(meancv)])
  }


gcvlam <- function(Gramat,y,mscale)
  { LL <- 10
    lam0 <- 2^(-(1:LL))
    Kmat <- wsGram(Gramat,mscale)
    gcvs <-  rep(0,LL)
    for (jj in 1:LL)
        {gcvs[jj] = gcv(Kmat,y,lam0[jj])}
    return(lam0[which.min(gcvs)])
  }
  
#input: x,y,lam,M
#output: cbtheta
twostep <- function(Gramat,y,wt,lam,mm)
  { n <- length(y)
    d <- dim(Gramat)[3]
    theta <- rep(1,d)
    cbtheta <- rep(0,n+d+1)
    cb0 <- sspline(Gramat,y,1/(wt^2),lam)$cb #calculate initial (c,b)
    c0 <- cb0[1:n]
    b0 <- cb0[1+n]
    G0 <- matrix(0,n,d)
    for (j in 1:d)
     {G0[,j] = Gramat[,,j]%*%c0*(wt[j]^(-2))}
    ycb0 = y-lam*c0/2-b0
    consM = matrix(0,d+1,d)
    consM[1:d,1:d] = diag(1,d)
    consM[d+1,] = -1 
    bvec = c(rep(0,d),-mm)
    #using constraint least squares (much slower than QP)
    #theta1 <- lsei(A=G0,B=ycb0,E=NULL,F=NULL,G=consM,H=bvec)[[1]]
    ## using QP.solve, recommended!
    Dmat <- 2*t(G0)%*%G0
    dvec <- 2*t(G0)%*%ycb0
    theta1 <- solve.QP(Dmat,dvec,t(consM),bvec)[[1]]
    sfit <- sspline(Gramat,y,theta1/(wt^2),lam)
    solution <- c(sfit$cb,theta1)
    Df<- sfit$Df 
    result <- list(solution=solution,Df=Df)
    return(result)
 }
  
  
  
## tune (lam,mm) using CV
cvadd <- function(Gramat,y,wt,folds)
{   n <- length(y)
    d <- dim(Gramat)[3]
    bestlam = cvlam(Gramat,y,wt^(-2),folds)
    splits = cvsplit(n,folds) 
    pp = sample(n,replace=FALSE) 
    aa <- seq(0.5,d*2,by=0.5)
    cverr <- matrix(0,length(aa),folds)
    for (kk in 1:folds)
    { tmp = 0*(kk==1)+sum(splits[1:(kk-1)])*(kk>1) #starting index 
      ntrain = n-splits[kk] #splits[kk] is the tuning size, ntrain = training size
      tuneid = pp[(tmp+1):(tmp+splits[kk])]
      ytune = y[tuneid]
      ytrain = y[-tuneid]
      trainGramat = Gramat[-tuneid,-tuneid,1:d]
      predGramat = Gramat[tuneid,-tuneid,1:d]
      for (jj in 1:length(aa))
        {jjcbtheta = twostep(trainGramat,ytrain,wt,bestlam,aa[jj])$solution
         jjcb = jjcbtheta[1:(ntrain+1)]
         jjtheta = jjcbtheta[(2+ntrain):(1+ntrain+d)]
         Kpred = wsGram(predGramat,jjtheta/(wt^2))
         jjpredict = Kpred%*%jjcb[1:ntrain]+jjcb[ntrain+1]
         cverr[jj,kk] = sum((ytune-jjpredict)^2)
        }
    }
      meancv = apply(cverr,1,sum)/n #the total leave-out error/n
      bestmm = aa[which.min(meancv)]
      result = twostep(Gramat,y,wt,bestlam,bestmm)
      return(result)
}

## bicvadd: tune lambda with GCV and tune mm with BIC (sigma2 is estimated by MLE)
bicadd <- function(Gramat,y,wt)
{   n <- length(y)
    d <- dim(Gramat)[3]
    bestlam = gcvlam(Gramat,y,1/(wt^2))
    if(d<30)   aa <- seq(0.5,d*2,by=0.5)
    else       aa <- seq(1,d,by=1)
    bics <- rep(0,length(aa))
    for (jj in 1:length(aa))
        {
         tmp = twostep(Gramat,y,wt,bestlam,aa[jj])
         jjcbtheta = tmp$solution
         jjdf = tmp$Df 
         jjc = jjcbtheta[1:n]
         jjb = jjcbtheta[1+n]
         jjtheta = jjcbtheta[(2+n):(1+n+d)]
         jjKmat = wsGram(Gramat,jjtheta/(wt^2))
         ypred =  jjKmat%*%jjc+jjb
         bics[jj] = bic(y,ypred,jjdf)
        }
      bestmm = aa[which.min(bics)]
      result = twostep(Gramat,y,wt,bestlam,bestmm)
      return(result)
}

##compute gcv for standard smoothing splines
gcv <- function(Kmat,y,lam)
{ n = length(y)
  myqr <- qr(matrix(1,n,1)) 
  myQ <- qr.Q(myqr,complete=TRUE)
  myR <- qr.R(myqr,complete=TRUE)
  F1 <- myQ[,1]
  F2 <- myQ[,-1]
  z <- t(F2)%*%y
  S <- t(F2)%*%Kmat%*%F2
  tmp = solve(S+lam*diag(1,n-1))
  tmpz = tmp%*%z
  trt = sum(diag(tmp))
  score = n*t(tmpz)%*%tmpz/(trt^2)
  return(score)
}


bic <- function(y,ypred,df)
{ n = length(y)
  rss = sum((y-ypred)^2)
  score = n*log(rss/n)+ log(n)*df
  return(score)
}


## compute df for smoothing splines
compdf <- function(Kmat,y,lam)
 {n = length(y)
  myqr <- qr(matrix(1,n,1)) 
  myQ <- qr.Q(myqr,complete=TRUE)
  myR <- qr.R(myqr,complete=TRUE)
  F1 <- myQ[,1]
  F2 <- myQ[,-1]
  S <- t(F2)%*%Kmat%*%F2
  tmp = solve(S+lam*diag(1,n-1))
  df = n-lam*sum(diag(tmp))
  return(df)
}

##compute weights
compwt <- function(x,y,mscale=rep(1,ncol(x)),gampow=1)
  { n <- length(y)
    d <- dim(as.matrix(x))[2]
    Gramat <- bigGram(x,x)
    bestlam = gcvlam(Gramat,y,mscale)  
    cc <- sspline(Gramat,y,mscale,bestlam)$cb[1:n]
    adwt <- rep(1,d)
    for (jj in 1:d)
     {fjj=Gramat[,,jj]%*%cc
      adwt[jj] <- (sqrt(mean(fjj^2)))^(-gampow)
     }
    return(adwt)
  }
