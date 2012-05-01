.packageName <- "cosso"
#######=============================#####
#---------- Common Functions -----------#
#######=============================#####

rescale <- function(x)
  {
  return( (x-min(x))/(max(x)-min(x)) )
  }

##generate spline kernel (including linear terms)
genK<- function(x1,x2)
    {
    m = length(x1)
    n = length(x2)
    b = abs(x1%*%matrix(1,ncol=n)-matrix(1,nrow=m)%*%t(x2))
    k1s = x1 - 0.5
    k1t = x2 - 0.5
    k2s = (k1s^2 - 1/12)/2
    k2t = (k1t^2 - 1/12)/2
    K =k1s%*%t(k1t)+k2s%*%t(k2t)-((b - 0.5)^4 - (b - 0.5)^2 / 2 + 7/240)/24
    }

gram <- function(x1,x2,mscale)
   { 
     n1 <- dim(as.matrix(x1))[1]
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
  { 
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  d <- ncol(x1)
  Gram <- array(0,c(n1,n2,d))
  for(j in 1:d)  Gram[,,j] <-  genK(x1[,j], x2[,j])
  return(Gram)
  }


wsGram <- function(Gramat,mscale)
   { 
    n1 <- dim(Gramat)[1]
    n2 <- dim(Gramat)[2]
    d <- dim(Gramat)[3]
    KK <- matrix(0,n1,n2)
    for (j in 1:d)   KK = KK + mscale[j]*Gramat[,,j]
    return(KK)
   }

predict.cosso<- function(object,xnew,M,type=c("fit","coefficients"),...)
  {
  type<- match.arg(type)
  if (missing(xnew) & type == "fit")
     {
      warning("Type=fit with no xnew argument; type switched to coefficients")
      type <- "coefficients"
     }
  if(missing(M))
    {
     warning("Missing Smoothing Parameter, M. Replaced by value tuned by 5-CV")
     if(is.null(object$tau))
        {  tuneobj<- tune.cosso(object,type="CV", folds=5,plot.it=FALSE)  }
     else
        {  tuneobj<- tune.cosso.qr(object,folds=5,plot.it=FALSE) }
     M <- tuneobj$OptM
    }
  n<- nrow(object$x)
  nbasis <- length(object$basis.id)
  d<- ncol(object$x)
  if(!is.null(object$tau))
     {
     solution  <- cquan(tau=object$tau,y=object$y,K3dtrain=object$Kmat,lam0=object$tune$OptLam,mm=M,wts=object$wt)
     }
  else
     {
     solutionset <- twostep(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,M)$solution
     solution <- list(intercept=solutionset[1+nbasis],coefs=solutionset[1:nbasis],theta=solutionset[-c(1:(nbasis+1))])
     }
  if (type=="fit")
     {
      Knew = gram(xnew, object$x[object$basis.id,], solution$theta/(object$wt^2))
      predictor = Knew%*%solution$coefs+solution$intercept
     }
  else
     {
      predictor = list(coef=solution$coefs,b=solution$intercept,theta=solution$theta)
     }
  class(predictor)="predict.cosso"
  return(predictor)
  }

plot.cosso<- function(x,M,plottype =c("Path","Functionals"),eps=1e-10,...)
  {
  object <- x
  d<- ncol(object$x)
  n<- nrow(object$x)
  nbasis<- length(object$basis.id)
  plottype <- match.arg(plottype)
  if (missing(M) & plottype == "Functionals")
     {
      warning("Type=fit with no smoothing parameter argument; type switched to Path")
      plottype <- "Path"
     }
  if(plottype=="Path")
     {
      matplot(object$tune$Mgrid,object$tune$L2norm,type="l",lty=1,col=c(1,rainbow(d-1)),xlab="M",ylab=expression(L[2]-norm))
      axis(4,at=object$tune$L2norm[length(object$tune$Mgrid),],labels=1:d,cex=.6,las=2)
     }
  else
     {
      if(!is.null(object$tau))
         {
         solution  <- cquan(tau=object$tau,y=object$y,K3dtrain=object$Kmat,lam0=object$tune$OptLam,mm=M,wts=object$wt)
         }
      else
         {
         solutionset <- twostep(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,M)$solution
         solution <- list(coefs=solutionset[1:nbasis],theta=solutionset[-c(1:(nbasis+1))])
         }
      sid <- (1:d)[abs(solution$theta)>eps]
      fsize <- length(sid)
      fits <- matrix(0,n,d)
      for (jj in 1:d)   fits[,jj]=solution$theta[jj]/(object$wt[jj]^2)*object$Kmat[,object$basis.id,jj]%*%solution$coefs
      frange = range(fits)
      plotrow <- ceiling(fsize/2)
      plotcol <- 2
      par(mfrow=c(plotrow,plotcol))
      for (jj in 1:fsize)  #only plot the selected components
        {
        sjj = sid[jj]
        jorder = order(object$x[,sjj])
        if (is.null(names(object$x))==0)
            plot(object$x[jorder,sjj],fits[jorder,sjj],type="l",xlim=c(0,1),ylim=frange,xlab=names(object$x)[sjj],ylab=paste("f(",names(object$x)[sjj],")",sep=""))
        else
          plot(object$x[jorder,sjj],fits[jorder,sjj],type="l",xlim=c(0,1),ylim=frange,xlab=paste("x",eval(sjj),sep=""),ylab=paste("f(x",eval(sjj),")",sep=""))
        }
     }
  }

#######=============================#####
#--- Functions for Gaussian Response ---#
#######=============================#####

cosso <- function(x,y,wt=rep(1,ncol(x)),scale=FALSE,nbasis,basis.id,n.step=2*ncol(x))
 { 
     n <- dim(as.matrix(x))[1]
     d <- dim(as.matrix(x))[2]
     if(scale) x <- apply(x,2,rescale)
     if( missing(nbasis) &  missing(basis.id))   { nbasis<- n ; basis.id <- 1:n }
     if( missing(nbasis) & !missing(basis.id))     nbasis <- n
     if(!missing(nbasis) &  missing(basis.id))     basis.id <- sort(sample(1:n,nbasis))
     GramatF <- bigGram(x,x)
     Gramat1 <- GramatF[,basis.id,]    
     Gramat2 <- GramatF[basis.id,basis.id,]    
     
     bestlam <- gcvlam(GramatF,y,1/(wt^2))
     L2normMat <- bicVec <- Mgrid <- NULL
     tempM<- 0.1
     tempTheta=rep(0,d)
     loop=0
     while(sum(tempTheta>1e-7)<d & loop<n.step)
         {
         loop=loop+1
         Mgrid <- c(Mgrid,tempM)
         tmp <- twostep(Gramat1,Gramat2,y,wt,bestlam,tempM)
         jjcbtheta <- tmp$solution
         jjc <- jjcbtheta[1:nbasis]
         jjb <- jjcbtheta[1+nbasis]
         tempTheta <- jjcbtheta[-c(1:(nbasis+1))]
         yhat <- jjb+wsGram(Gramat1,tempTheta/wt^2)%*%jjc
         bicVec <- c(bicVec,bic(y,yhat,tmp$Df))
         #--- Compute Solution Path ---#
         for(j in 1:d)  L2normMat<- c( L2normMat,sqrt(mean((tempTheta[j]/wt[j]^2*Gramat1[,,j]%*%jjc)^2)) )
         tempM=tempM+0.5
         }
     L2normMat=matrix(L2normMat,ncol=d,byrow=TRUE)
     L2normMat=rbind(rep(0,d),L2normMat)
     Mgrid=c(0,Mgrid)
     bicVec=c(NA,bicVec)
     cossoobj<- list( x=x,y=y,wt=wt,Kmat=GramatF,basis.id=basis.id,tune=list(OptLam=bestlam,bic=bicVec,Mgrid=Mgrid,L2norm=L2normMat) )
     class(cossoobj)="cosso"
     return(cossoobj)
 }



tune.cosso <- function(object,type=c("BIC","CV"), folds=5,plot.it=TRUE)
   {
     n <- length(object$y)
     nbasis <- length(object$basis.id)
     d <- dim(as.matrix(object$x))[2]
     type <- match.arg(type)
     if (type=="BIC")
       {tuneobj <- bicadd(object)}
     if (type=="CV")
       {tuneobj <- cvadd(object,folds)}
     if(plot.it)
       {
       par(mfcol=c(1,2))
       plotid=complete.cases(cbind(tuneobj$Mgrid,tuneobj$IC))
       plot(tuneobj$Mgrid[plotid],tuneobj$IC[plotid],type="l",lwd=1.5,xlab="M",ylab=type)
       abline(v=tuneobj$OptM,lty=2,col=2);axis(3,tuneobj$OptM)
       matplot(tuneobj$Mgrid,tuneobj$L2norm,type="l",lty=1,col=c(1,rainbow(d-1)),xlab="M",ylab=expression(L[2]-norm))
       abline(v=tuneobj$OptM,lty=2,col=2);axis(3,tuneobj$OptM)
       axis(4,at=tuneobj$L2norm[length(tuneobj$Mgrid),],labels=1:d,cex=.3,las=2)
       }
     return(tuneobj)
    }

## solve the standard smoothing spline
sspline <- function(Gramat1,Gramat2,y,mscale,lam)
  { 
    n <- length(y)
    d <- length(mscale)
    nbasis <- dim(Gramat1)[2]
    if(nbasis==n) 
      {
      Kmat = wsGram(Gramat1,mscale)
      bigK = matrix(1,n+1,n+1)
      bigK[1:n,1:n] = Kmat+lam*diag(1,n)
      bigK[n+1,n+1]= 0
      cb = solve(bigK,c(y,0))
      Df = compdf(Kmat,y,lam)
      }
    else
      {
      Rtheta1 = wsGram(Gramat1,mscale)
      Rtheta2 = wsGram(Gramat2,mscale)
      Mprime=solve(t(Rtheta1)%*%Rtheta1+lam*Rtheta2+1e-8*diag(dim(Rtheta2)[1]))%*%t(Rtheta1)
      M=Rtheta1%*%Mprime
      c1=sum(diag(n)-M)
      bhat=sum((diag(n)-M)%*%y)/c1
      chat=Mprime%*%(diag(n)-matrix(1,ncol=n,nrow=n)%*%(diag(n)-M)/c1)%*%y
      cb=c(chat,bhat)
      Smat=M+(diag(n)-M)%*%matrix(1,ncol=n,nrow=n)%*%(diag(n)-M)/c1
      Df=sum(diag(Smat))    
      }
    result=list(cb=cb,Df=Df)
    return(result)
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

twostep <- function(Gramat1,Gramat2,y,wt,lam,mm)
   { 
    n <- length(y)
    nbasis <- dim(Gramat1)[2]
    d <- length(wt)

    cb0<- sspline(Gramat1,Gramat2,y,1/(wt^2),lam)$cb
    c0 <- cb0[1:nbasis]
    b0 <- cb0[1+nbasis]
    #-------------------------------------#
    G1 <- matrix(0,n,d)
    G2 <- matrix(0,nbasis,d)
    for(j in 1:d)
       {
       G1[,j] = Gramat1[,,j]%*%c0*(wt[j]^(-2))
       G2[,j] = Gramat2[,,j]%*%c0*(wt[j]^(-2))
       }
    dvec <- 2*t(G1)%*%(y-b0)-lam*t(G2)%*%c0
    Dmat <- 2*t(G1)%*%G1
    Amat <- rbind(diag(d),rep(-1,d))
    bvec   <- c(rep(0,d),-mm)
    theta1 <- solve.QP(Dmat,dvec,t(Amat),bvec)[[1]]
    theta1[theta1<1e-8] <- 0
    #-------------------------------------#
    cb1 <- sspline(Gramat1,Gramat2,y,theta1/(wt^2),lam)
    solution <- c(cb1$cb,theta1)
    result <- list(solution=solution,Df=cb1$Df)
    return(result)
  }


cvadd <- function(object,folds)
   {
    n <- length(y)
    d <- dim(Gramat)[3]
    Gramat <- object$Kmat
    y <- object$y
    wt<- object$wt
    if(d<20)   mm <- seq(0.25,ceiling(max(object$tune$Mgrid)),by=0.5)
    else       mm <- seq(0.25,ceiling(max(object$tune$Mgrid)),by=1)
    Fullmm=as.numeric(names(table(c(mm,object$tune$Mgrid))))
    
    splits =  cvsplit(n,folds)
    pp = sample(n,replace=FALSE)
    cverr <- matrix(0,length(Fullmm),folds)
    cvVec <- rep(NA,length(Fullmm))
    L2normMat <- matrix(NA,ncol=d,nrow=length(Fullmm))
    prevTheta=rep(0,d)
    for (jj in 1:length(Fullmm))
      {
      if(!Fullmm[jj]%in%object$tune$Mgrid)
          {
          #--- Without CV. Compute Solution Path ---#
          coefhat   <- twostep(Gramat,Gramat,y,wt,object$tune$OptLam,Fullmm[jj])$solution
          tempc     <- coefhat[1:n]
          temptheta <- coefhat[(2+n):(1+n+d)]
          for(j in 1:d)  L2normMat[jj,j]=sqrt(mean((temptheta[j]/wt[j]^2*Gramat[,,j]%*%tempc)^2))
    
          #----- CV ------#
          for (kk in 1:folds)
              {
              tmp = 0*(kk==1)+sum(splits[1:(kk-1)])*(kk>1) #starting index
              ntrain = n-splits[kk] #splits[kk] is the tuning size, ntrain = training size
              tuneid = pp[(tmp+1):(tmp+splits[kk])]
              ytune = y[tuneid]
              ytrain = y[-tuneid]
              trainGramat = Gramat[-tuneid,-tuneid,]
              predGramat  = Gramat[tuneid,-tuneid,]
      
              jjcbtheta = twostep(trainGramat,trainGramat,ytrain,wt,object$tune$OptLam,Fullmm[jj])$solution
              jjcb = jjcbtheta[1:(ntrain+1)]
              jjtheta = jjcbtheta[(2+ntrain):(1+ntrain+d)]
              Kpred = wsGram(predGramat,jjtheta/(wt^2))
              jjpredict = Kpred%*%jjcb[1:ntrain]+jjcb[ntrain+1]
              cverr[jj,kk] = sum((ytune-jjpredict)^2)
              }
           cvVec[jj]=mean(cverr[jj,])
           }   
      }
    bestmm = Fullmm[which.min(cvVec)]
    L2normMat[is.na(L2normMat[,1]),]=object$tune$L2norm
    tuneobj=list(OptM=bestmm,OptLam=object$tune$OptLam,Mgrid=Fullmm,IC=cvVec,L2norm=L2normMat)
    return(tuneobj)
   }


bicadd=function(object)
  {
  n=length(object$y)
  nbasis=length(object$basis.id)
  d=ncol(object$x)
  origMgrid=object$tune$Mgrid
  refinePt=which(apply(object$tune$L2norm<1e-8,1,sum)[-length(origMgrid)]-apply(object$tune$L2norm<1e-8,1,sum)[-1]>1)
  if(length(refinePt)==0)
     {   extMgrid=c(apply(cbind(origMgrid[-1],origMgrid[-length(origMgrid)]),1,mean),max(origMgrid)+0.5)  }
  else
     {   extMgrid=c(as.numeric(apply(cbind(origMgrid[refinePt],origMgrid[refinePt+1]),1,quantile,c(.2,.4,.6,.8))),max(origMgrid)+0.5,d)  }
  extbicVec=rep(NA,length(extMgrid))
  extL2normMat=matrix(NA,ncol=d,nrow=length(extMgrid))
  for(jj in 1:length(extMgrid))
     {
      tmp = twostep(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,extMgrid[jj])
      jjcbtheta = tmp$solution
      jjdf = tmp$Df
      jjc = jjcbtheta[1:nbasis]
      jjb = jjcbtheta[1+nbasis]
      jjtheta = jjcbtheta[-c(1:(nbasis+1))]
      #--- Compute Solution Path ---#
      for(j in 1:d)  extL2normMat[jj,j]=sqrt(mean((jjtheta[j]/object$wt[j]^2*object$Kmat[,object$basis.id,j]%*%jjc)^2))
      #--- Compute BIC ---#
      jjKmat = wsGram(object$Kmat[,object$basis.id,],jjtheta/(object$wt^2))
      ypred =  jjKmat%*%jjc+jjb
      extbicVec[jj] = bic(object$y,ypred,jjdf)
      }
   Mgrid=c(origMgrid,extMgrid)
   bicVec=c(object$tune$bic,extbicVec)[order(Mgrid)]
   L2norm=rbind(object$tune$L2norm,extL2normMat)[order(Mgrid),]
   Mgrid=Mgrid[order( Mgrid)]
   bestmm=Mgrid[which.min(bicVec)]
   tuneobj=list(OptM=bestmm,OptLam=object$tune$OptLam,Mgrid=Mgrid,IC=bicVec,L2norm=L2norm)
   return(tuneobj)
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
SSANOVAwt <- function(x,y,mscale=rep(1,ncol(x)),gampow=1)
  { n <- length(y)
    d <- dim(as.matrix(x))[2]
    Gramat <- bigGram(x,x)
    bestlam = gcvlam(Gramat,y,mscale)
    cc <- sspline(Gramat,Gramat,y,mscale,bestlam)$cb[1:n]
    adwt <- rep(1,d)
    for (jj in 1:d)
     {fjj=Gramat[,,jj]%*%cc
      adwt[jj] <- (sqrt(mean(fjj^2)))^(-gampow)
     }
    return(adwt)
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

#######=============================#####
#-- Functions for Quantile Regression --#
#######=============================#####

cosso.qr=function(x,y,tau,wt=rep(1,ncol(x)),scale=FALSE,parallel=FALSE,cpus=1)
    {
     n<- nrow(x)
     d<- ncol(x)
     if(scale) x <- apply(x,2,rescale)
     Gramat  <- bigGram(x,x)
     #--- Initializate Parallel Computing ---#
     sfInit(parallel=parallel,cpus=cpus)
     bestlam <- cvLam.qr(Gramat,y,tau,wt,folds=5,parallel=parallel,cpus=cpus)
     tempM <- seq(0.5,d*0.6,0.75)
     L2normMat <- matrix(NA,ncol=d,nrow=length(tempM))

     if(parallel) 
        { 
         sfLibrary(cosso);sfLibrary(quadprog);sfLibrary(Rglpk)  
        }
     tempcoefs=sfClusterApplyLB(tempM,cquan,tau=tau,y=y,K3dtrain=Gramat,lam0=bestlam,wts=wt)
     for(m in 1:length(tempM))
        {
        coefhat=tempcoefs[[m]]
        for(j in 1:d)  L2normMat[m,j]=sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat[,,j]%*%coefhat$coefs)^2))
        }
        
     if(sum(L2normMat[length(tempM),]==0)>0)  # Some component remain unselected
       {
       extMgrid=seq(max(tempM)+1,d*0.85,l=5)
       extL2normMat=matrix(NA,ncol=d,nrow=length(extMgrid))
       tempcoefs=sfClusterApplyLB(extMgrid,cquan,tau=tau,y=y,K3dtrain=Gramat,lam0=bestlam,wts=wt)
       for(m in 1:length(extMgrid))
          {
          coefhat=tempcoefs[[m]]
          for(j in 1:d)  extL2normMat[m,j]=sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat[,,j]%*%coefhat$coefs)^2))
          }
       tempM=c(tempM,extMgrid)   
       L2normMat=rbind(L2normMat,extL2normMat)
       }
     sfStop()

     cossoqrobj<- list(x=x,y=y,tau=tau,wt=wt,Kmat=Gramat,basis.id=1:n,tune=list(OptLam=bestlam,Mgrid=c(0,tempM),L2norm=rbind(rep(0,d),L2normMat)) )
     class(cossoqrobj)="cosso"
     return(cossoqrobj)
    }

tune.cosso.qr <- function(object,folds=5,plot.it=TRUE,parallel=FALSE,cpus=1)
   {
     tuneobj <- CVadd(object,folds,parallel,cpus)
     if(plot.it)
       {
       d<- length(object$wt)
       validID<- complete.cases(tuneobj$Mgrid,tuneobj$L2norm[,1])
       par(mfcol=c(1,2))
       plot(tuneobj$Mgrid,tuneobj$CV,type="l",lwd=1.5,xlab="M",ylab="CV score")
       abline(v=tuneobj$OptM,lty=2,col=2);axis(3,tuneobj$OptM)
       matplot(tuneobj$Mgrid[validID],tuneobj$L2norm[validID,],type="l",lty=1,col=c(1,rainbow(d-1)),xlab="M",ylab=expression(L[2]-norm))
       abline(v=tuneobj$OptM,lty=2,col=2);axis(3,tuneobj$OptM)
       axis(4,at=tuneobj$L2norm[length(tuneobj$Mgrid),],labels=1:d,cex=.3,las=2)
       }
     return(tuneobj)
    }

cquan=function(tau,y,K3dtrain,lam0,mm,wts)
  {
  n <- length(y)
  d <- length(wts)
  G <- matrix(0,n,d)

  newtheta<- rep(1,d)     
  Ktheta0 <- wsGram(K3dtrain,newtheta/wts^2)
  cb   <- kqr(y,tau,lam0,Ktheta0,insure=TRUE)
  newc <- cb$coefs
  newb <- as.numeric(cb$intercept)

  for (i in 1:d)  G[,i]=(wts[i]^(-2))*K3dtrain[,,i]%*%as.vector(newc)
  newtheta <- garrote.qr(x=G,y=y-newb,tau=tau,lambda0=lam0,M=mm,ct=t(newc))

  Ktheta0 <- wsGram(K3dtrain,newtheta/wts^2)
  cb <- kqr(y,tau,lam0,Ktheta0,insure=TRUE)
  newc <- cb$coefs
  newb <- as.numeric(cb$intercept)

  output<-list(coefs=newc,intercept=newb,theta=newtheta,quantile=cb$quantile)
  return(output)
  }


kqr=function(y,tau,lambda,Ktheta,insure=TRUE,posconst=1e-8)
   {
    n<-length(y)
    r<-1/(2*n*lambda)
    K<-Ktheta
    if(insure)      K<-K+posconst*diag(1,n)
    #to guarantee a positive definite kernel matrix
    Amat<-rbind(diag(rep(1,n)),diag(rep(-1,n)))
    Amat<-t(rbind(rep(1,n),Amat))
    b0<-c(0,rep(r*(tau-1),n),rep(-r*tau,n))
    hatalpha<-drop(solve.QP(Dmat=K,dvec=y,Amat=Amat,bvec=b0,meq=1)$solution)
    hatf<-drop(K%*%hatalpha)
    hatbeta0<-quantile(y-hatf,tau)
    output<-list(coefs=hatalpha,intercept=hatbeta0,quantile=(hatf+hatbeta0))
    output
   }


garrote.qr=function(x,y,tau,lambda0,ct,M)
   {
   #minimize 1/n sum_i rho_tau(y_i-G_i theta)+lambda0 c^t G theta  subject to sum_j theta_j <=M, theta_j>=0
   #x is the G matrix
   n<-length(y)
   p<-dim(x)[2]
   cvec<-c(lambda0*drop(ct%*%x)-(tau-0.5)*apply(x,2,mean),rep(1,n)/(2*n))
   mat.con<-cbind(diag(p),matrix(0,p,n))
   mat.con<-rbind(mat.con,c(-rep(1,p),rep(0,n)))
   mat.con<-rbind(mat.con,cbind(x,diag(n)))
   mat.con<-rbind(mat.con,cbind(-x,diag(n)))
   lp.out<-Rglpk_solve_LP(obj=cvec,mat=mat.con,dir=rep(">=",2*n+p+1),rhs=c(rep(0,p),-M,y,-y),max=FALSE)$solution
   theta<-lp.out[1:p]
   theta[theta<1e-5]<-0
   return(theta)
   }

CVadd=function(object,folds,parallel,cpus)
    {
    d<- length(object$wt)
    #--  Tuning Grids --$
    origMgrid<- object$tune$Mgrid
    refinePt<- which(apply(object$tune$L2norm<1e-8,1,sum)[-length(origMgrid)]-apply(object$tune$L2norm<1e-8,1,sum)[-1]>1)
    if(length(refinePt)>0)
       {   extMgrid<- c( as.numeric(apply(cbind(origMgrid[refinePt],origMgrid[refinePt+1]),1,mean)) )  }
    else
       {   extMgrid<- NULL }
    candmm <- c(origMgrid,extMgrid)   
    cvMobj <- cvM.qr(object$Kmat,object$y,object$tau,object$wt,candmm[-1],object$tune$OptLam,folds,parallel,cpus)
    
    cvscore<- c(NA,cvMobj$CVerr)[order(candmm)]
    extL2norm<- matrix(NA,ncol=d,nrow=length(extMgrid))
    L2norm<- rbind(object$tune$L2norm,extL2norm)[order(candmm),]
    candmm<- sort(candmm)

    tuneobj=list(OptM=cvMobj$OptM,OptLam=object$tune$OptLam,Mgrid=candmm,CV=cvscore,L2norm=L2norm)
    return(tuneobj)
    }

cvLam.qr=function(Gramat,y,tau,wt,lam,folds,theta,parallel=FALSE,cpus=1)
   {
   n<- length(y); d<- length(wt)
   if(missing(theta)) theta <- rep(1,d)
   if(missing(lam))   lam   <- 2^seq(-15,-23,by=-0.5)

   CVmat  <- matrix(NA,ncol=length(lam),nrow=folds)
   splitID<- cvsplitID(n,folds)

   sfLibrary(quadprog)
   for(f in 1:folds)
      {
        testID <- splitID[!is.na(splitID[,f]),f]
        trainID<- (1:n)[-testID]
        trainRtheta<- wsGram(Gramat[trainID,trainID,],theta/wt^2)
        testRtheta <- wsGram(Gramat[testID ,trainID,],theta/wt^2)

        coefhat<- sfClusterApplyLB(lam,kqr,Ktheta=trainRtheta,y=y[trainID],tau=tau)

        for(l in 1:length(lam))
           {
           yhat   <- testRtheta%*%coefhat[[l]]$coefs+coefhat[[l]]$intercept
           CVmat[f,l]<- sum(rho(tau, y[testID]-yhat))
           }
      }

   optLam<- lam[which.min(apply(CVmat,2,sum))]
   return(optLam)
   }


 cvM.qr<- function(Gramat,y,tau,wt,M,lam0,folds,parallel=FALSE,cpus=1)
   {
   n<- length(y); d<- length(wt)
   CVmat  <- matrix(NA,ncol=length(M),nrow=folds)
   splitID<- cvsplitID(n,folds)

   sfInit(parallel=parallel, cpus=cpus)
   sfLibrary(quadprog);sfLibrary(Rglpk)
   if(parallel) sfLibrary(cosso)
   for(f in 1:folds)
      {
        testID <- splitID[!is.na(splitID[,f]),f]
        trainID<- (1:n)[-testID]
        trainK3darray<- Gramat [trainID,trainID,]
        testK3darray <- Gramat [testID ,trainID,]
        #--- Parallel ---#
        coefhat<- sfClusterApplyLB(M,cquan,tau=tau,y=y[trainID],K3dtrain=trainK3darray,lam0=lam0,wts=wt)

        for(mm in 1:length(M))
           {
            testRtheta<- wsGram(testK3darray,coefhat[[mm]]$theta/wt^2)
            yhat<- testRtheta%*%coefhat[[mm]]$coefs+coefhat[[mm]]$intercept
            CVmat[f,mm]<- sum(rho(tau, y[testID]-yhat))
           }
      }
   sfStop()
   optM <- M[which.min(apply(CVmat,2,sum))]
   return(list(OptM=optM,CVerr=apply(CVmat,2,sum)/n))
   }

 
cvsplitID=function(n,folds)
  {
    fsize  <- floor(n/folds) #average size of each fold, the last fold is larger
    splits <- fsize*rep(1,folds)
    nextra <- n-folds*fsize
    if (nextra>0)
      {splits[1:nextra] <- splits[1:nextra]+1}
    randid<- sample(1:n,n)
    IDmat <- matrix(NA,ncol=folds,nrow=ceiling(n/folds))
    IDmat[,1]<- randid[1:splits[1]]
    for(i in 2:folds)
       {
        tempid<-randid[(cumsum(splits)[i-1]+1):(cumsum(splits)[i])]
        length(tempid)<- ceiling(n/folds)
        IDmat[,i]<- tempid
       }
  return(IDmat)
  }

df.est=function(object,tol=1e-7)
   {
   dfhat<- sum(abs(object$y-object$quantile)<tol)
   return(dfhat)
   }


rho<- function(tau,r)
   {
   sol <- tau*r*(r>0)-(1-tau)*r*(r<=0)
   return(sol)
   }

KQR=function(x,y,tau,lam0s,folds=5,parallel=FALSE,cpus=1)
    {
    n<- length(y); d<- ncol(x)

    Gram   <- bigGram(x,x)
    Rtheta <- wsGram(Gram,rep(1,d))

    if(missing(lam0s))  lam0s<- 2^seq(-14,-20,by=-0.5)
    TuneErrMat <- matrix(NA,nrow=folds,ncol=length(lam0s))
    splitID    <- cvsplitID(n,folds)
    
    sfInit(parallel=parallel, cpus=cpus)
    sfLibrary(quadprog)
    for(f in 1:folds)
        {
        testID <- splitID[!is.na(splitID[,f]),f]
        trainID<- (1:n)[-testID]
        trainRtheta<- wsGram(Gram[trainID,trainID,],rep(1,d))
        testRtheta <- wsGram(Gram[testID ,trainID,],rep(1,d))
        #--- Parallel Computing ---#
        coefhat<- sfClusterApplyLB(lam0s,kqr,Ktheta=trainRtheta,y=y[trainID],tau=tau)

        for(l in 1:length(lam0s))
           {
           yhat   <- as.numeric(testRtheta%*%coefhat[[l]]$coefs+coefhat[[l]]$intercept)
           TuneErrMat[f,l]<- sum(rho(tau, y[testID]-yhat))
           }
       }
     sfStop()
     optLam=lam0s[which.min(apply(TuneErrMat,2,sum))]

    result=kqr(y=y,tau=tau,lambda=optLam,Ktheta=Rtheta,insure=TRUE)
    
    kqrnorms=rep(NA,d)
    for(pp in 1:d)    kqrnorms[pp]=sqrt(mean((Gram[,,pp]%*%result$coef)^2))
    
    return(c(result,list(L2norm=kqrnorms,OptLam=optLam)))
    }


KQRwt<- function(x,y,tau,gampow=1,folds=5,parallel=FALSE,cpus=1)
  {
  KQRfit=KQR(x,y,tau,folds=folds,parallel=parallel,cpus=cpus)
  estiwt=1/KQRfit$L2norm^gampow
  return(estiwt)
  }
