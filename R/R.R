.packageName <- "cosso"
#######=============================#####
#---------- Common Functions -----------#
#######=============================#####

rescale <- function(x)
  {
  return( (x-min(x))/(max(x)-min(x)) )
  }

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

genK.cat <- function(x1,x2)
  {
  n1 <- length(x1)
  n2 <- length(x2)
  x1 <- rep(x1, times=n2)
  x2 <- rep(x2, each =n1)
  L  <- length(unique(x1))
  K  <- matrix(L*(x1==x2) - 1, n1, n2)
  return(K)
  }

bigGram <- function(x1,x2)
  {
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  d <- ncol(x1)
  Gram <- array(0,c(n1,n2,d))
  for(j in 1:d)
    {
    if(length(unique(x1[,j]))>6)   Gram[,,j] <-  genK(x1[,j], x2[,j])
    else                           Gram[,,j] <-  genK.cat(x1[,j], x2[,j])
    }
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


#==================================#
#------- Generic Functions --------#
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
     if(object$family=="Gaussian")
        {  tuneobj<- tune.cosso(object,type="CV", folds=5,plot.it=FALSE)  }
     else if(object$family=="Quantile")
        {  tuneobj<- tune.cosso.qr(object,folds=5,plot.it=FALSE) }
     else if(object$family=="Cox")   
        {  tuneobj<- tune.cosso.cox(object,plot.it=FALSE) }
     M <- tuneobj$OptM
    }
  n<- nrow(object$x)
  nbasis <- length(object$basis.id)
  d<- ncol(object$x)
  if(object$family=="Gaussian")
     {
     solutionset <- twostep(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,M)$solution
     solution <- list(intercept=solutionset[1+nbasis],coefs=solutionset[1:nbasis],theta=solutionset[-c(1:(nbasis+1))])
     }
  else if(object$family=="Quantile")
     {
     solution  <- cquan(tau=object$tau,y=object$y,K3dtrain=object$Kmat,lam0=object$tune$OptLam,mm=M,wts=object$wt)
     }
  else if(object$family=="Cox") 
     {
     
     }
  if (type=="fit")
     {
      Knew = wsGram(bigGram(xnew, object$x[object$basis.id,]), solution$theta/(object$wt^2))
      predictor = Knew%*%solution$coefs+solution$intercept
     }
  else
     {
      predictor = list(coef=solution$coefs,b=solution$intercept,theta=solution$theta)
     }
  class(predictor)="predict.cosso"
  return(predictor)
  }

plot.cosso<- function(x,M,plottype =c("Path","Functionals"),eps=1e-7,...)
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
      if(object$family=="Gaussian")
         {
         solutionset <- twostep(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,M)$solution
         solution <- list(coefs=solutionset[1:nbasis],theta=solutionset[-c(1:(nbasis+1))])
         }
      else if(object$family=="Quantile")
         {
         solution  <- cquan(tau=object$tau,y=object$y,K3dtrain=object$Kmat,lam0=object$tune$OptLam,mm=M,wts=object$wt)
         }
      else if(object$family=="Cox")
         {
         solution <- twostep.Cox(Gramat1=object$Kmat[,object$basis.id,],Gramat2=object$Kmat[object$basis.id,object$basis.id,],time=object$time,status=object$status,wt=object$wt,basis.id=object$basis.id,RS=object$RiskSet,lambda0=object$tune$OptLam,M=M)
         }
      sid <- (1:d)[abs(solution$theta)>eps]
      fsize <- length(sid)
      fits <- matrix(0,n,d)
      for (jj in 1:d)   fits[,jj]=solution$theta[jj]/(object$wt[jj]^2)*object$Kmat[,object$basis.id,jj]%*%solution$coefs
      plotrow <- ceiling(fsize/2)
      plotcol <- 2
      par(mfrow=c(plotrow,plotcol))
      for (jj in 1:fsize)  #only plot the selected components
        {
        sjj = sid[jj]
        jorder = order(object$x[,sjj])
        if (is.null(names(object$x))==0)
            plot(object$x[jorder,sjj],fits[jorder,sjj],type="l",xlim=c(0,1),xlab=names(object$x)[sjj],ylab=paste("f(",names(object$x)[sjj],")",sep=""))
        else
          plot(object$x[jorder,sjj],fits[jorder,sjj],type="l",xlim=c(0,1),xlab=paste("x",eval(sjj),sep=""),ylab=paste("f(x",eval(sjj),")",sep=""))
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
     if( missing(nbasis) & !missing(basis.id))     nbasis <- length(basis.id)
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
     cossoobj<- list(family="Gaussian",x=x,y=y,wt=wt,Kmat=GramatF,basis.id=basis.id,tune=list(OptLam=bestlam,bic=bicVec,Mgrid=Mgrid,L2norm=L2normMat) )
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

    cb1 <- sspline(Gramat1,Gramat2,y,theta1/(wt^2),lam)
    solution <- c(cb1$cb,theta1)
    result <- list(solution=solution,Df=cb1$Df)
    return(result)
  }


cvadd <- function(object,folds)
   {
    y <- object$y
    Gramat <- object$Kmat
    wt<- object$wt
    n <- length(y)
    d <- dim(Gramat)[3]
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
     {   extMgrid <- c(apply(cbind(origMgrid[-1],origMgrid[-length(origMgrid)]),1,mean)[2:ceiling(length(origMgrid)/3+1)],max(origMgrid)+0.5)    }
  else
     {   extMgrid <- c(as.numeric(apply(cbind(origMgrid[refinePt],origMgrid[refinePt+1]),1,quantile,c(.3,.6)))          ,max(origMgrid)+0.5)    }
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
         sfLibrary("cosso",character.only=TRUE)
         sfLibrary("quadprog",character.only=TRUE)
         sfLibrary("Rglpk",character.only=TRUE)  
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

     cossoqrobj<- list(family="Quantile",x=x,y=y,tau=tau,wt=wt,Kmat=Gramat,basis.id=1:n,tune=list(OptLam=bestlam,Mgrid=c(0,tempM),L2norm=rbind(rep(0,d),L2normMat)) )
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

   if(parallel) sfLibrary("quadprog",character.only=TRUE)
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
   if(parallel) 
      {
      sfLibrary("cosso",character.only=TRUE)
      sfLibrary("quadprog",character.only=TRUE)
      sfLibrary("Rglpk",character.only=TRUE)
      }

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


KQRwt<- function(x,y,tau,cand.lam0,folds=5,parallel=FALSE,cpus=1)
  {
    n<- length(y); d<- ncol(x)

    Gram   <- bigGram(x,x)
    Rtheta <- wsGram(Gram,rep(1,d))

    if(missing(cand.lam0))  cand.lam0<- 2^seq(-14,-20,by=-0.75)
    TuneErrMat <- matrix(NA,nrow=folds,ncol=length(cand.lam0))
    splitID    <- cvsplitID(n,folds)
    
    sfInit(parallel=parallel, cpus=cpus)
    sfLibrary("quadprog",character.only=TRUE)
    for(f in 1:folds)
        {
        testID <- splitID[!is.na(splitID[,f]),f]
        trainID<- (1:n)[-testID]
        trainRtheta<- wsGram(Gram[trainID,trainID,],rep(1,d))
        testRtheta <- wsGram(Gram[testID ,trainID,],rep(1,d))
        #--- Parallel Computing ---#
        coefhat<- sfClusterApplyLB(cand.lam0,kqr,Ktheta=trainRtheta,y=y[trainID],tau=tau)

        for(l in 1:length(cand.lam0))
           {
           yhat   <- as.numeric(testRtheta%*%coefhat[[l]]$coefs+coefhat[[l]]$intercept)
           TuneErrMat[f,l]<- sum(rho(tau, y[testID]-yhat))
           }
       }
    sfStop()
    optLam=cand.lam0[which.min(apply(TuneErrMat,2,sum))]

    result=kqr(y=y,tau=tau,lambda=optLam,Ktheta=Rtheta,insure=TRUE)
    
    kqrnorms=rep(NA,d)
    for(pp in 1:d)    kqrnorms[pp]=sqrt(mean((Gram[,,pp]%*%result$coef)^2))
    return(1/kqrnorms)
  }

#######=============================#####
#------ Functions for Cox PH Model -----#
#######=============================#####
My_solve.QP=function(Dmat,dvec,Amat,bvec)
  {
  solution=tryCatch(solve.QP(Dmat,dvec,Amat,bvec)$solution,error=function(x) NA)
  if(is.na(solution[1]))
     {
     Dmat=diag(diag(Dmat))
     solution=solve.QP(Dmat,dvec,Amat,bvec)$solution
     }
   return(solution)
  }
My_solve=function(A,b)
  {
  solution=tryCatch(solve(A,b),error=function(x) NA)
  if(is.na(solution[1]))
     {
     solution=b/diag(A)
     }
   return(solution)
  }


cosso.cox <- function(x,time,status,wt=rep(1,ncol(x)),scale=FALSE,nbasis,basis.id,parallel=FALSE,cpus=1)
    {
     n<- nrow(x)
     d<- ncol(x)
     if(scale) x <- apply(x,2,rescale)

     if(missing(nbasis) & missing(basis.id))
       {
       nbasis=min( max(35, ceiling(12 * length(time)^(2/9)))  ,  sum(status==1)-5  )
       basis.id=sort(sample(which(status==1),nbasis))
       }
     else if(!missing(nbasis) & missing(basis.id))
       {
       if(nbasis>sum(status==1)) { nbasis=sum(status==1)-5; cat("Use nbasis=",sum(status==1)-5,"\n")}
       basis.id=sort(sample(which(status==1),nbasis))
       }

     Gramat <-bigGram(x,x)
     Gramat1<-Gramat[,basis.id, ]
     Gramat2<-Gramat1[basis.id,,]
     RS<-RiskSet(time,status)

     bestlam <- ACV.lambda(Gramat,time,status,1/wt^2,basis.id,RS)
     #--- Initializate Parallel Computing ---#
     sfInit(parallel=parallel,cpus=cpus)
     if(d<=10)  tempM <- seq(0.5,d*0.7,0.6)
     else       tempM <- c(seq(1,8,.75),seq(9,d*.7,1))
     L2normMat <- matrix(NA,ncol=d,nrow=length(tempM))
     ACVscore  <- rep(NA,length(tempM))
     if(parallel)
        {
        sfLibrary("quadprog",character.only=TRUE);sfLibrary("glmnet",character.only=TRUE)
        sfLibrary("cosso",character.only=TRUE)
        }
     tempcoefs <- sfClusterApplyLB(tempM,twostep.Cox,Gramat1=Gramat1,Gramat2=Gramat2,time=time,status=status,wt=wt,basis.id=basis.id,RS=RS,lambda0=bestlam)

     for(m in 1:length(tempM))
        {
        coefhat <- tempcoefs[[m]]
        ACVscore[m] <- PartialLik(time,status,RS,coefhat$fit)+sum(status==1)/n^2*( sum(diag(coefhat$UHU))/(n-1) - sum(coefhat$UHU)/(n^2-n) )
        for(j in 1:d)  L2normMat[m,j] <- sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat[,basis.id,j]%*%coefhat$coefs)^2))
        }

     if(sum(L2normMat[length(tempM),]==0)>0)  # Some component remain unselected
       {
       extMgrid=seq(max(tempM)+0.5,d*0.85,l=5)
       extL2normMat <- matrix(NA,ncol=d,nrow=length(extMgrid))
       extACVscore  <- rep(NA,length(extMgrid))
       tempcoefs=sfClusterApplyLB(extMgrid,twostep.Cox,Gramat1=Gramat1,Gramat2=Gramat2,time=time,status=status,wt=wt,basis.id=basis.id,RS=RS,lambda0=bestlam)
       for(m in 1:length(extMgrid))
          {
          coefhat=tempcoefs[[m]]
          extACVscore[m] <- PartialLik(time,status,RS,coefhat$fit)+sum(status==1)/n^2*( sum(diag(coefhat$UHU))/(n-1) - sum(coefhat$UHU)/(n^2-n) )
          for(j in 1:d)  extL2normMat[m,j] <- sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat[,basis.id,j]%*%coefhat$coefs)^2))
          }
       tempM=c(tempM,extMgrid)
       ACVscore <- c(ACVscore,extACVscore)
       L2normMat=rbind(L2normMat,extL2normMat)
       }
     sfStop()

     cossoqrobj<- list(family="Cox",x=x,time=time,status=status,wt=wt,Kmat=Gramat,basis.id=basis.id,RiskSet=RS,tune=list(OptLam=bestlam,ACV=c(NA,ACVscore),Mgrid=c(0,tempM),L2norm=rbind(rep(0,d),L2normMat)) )
     class(cossoqrobj)="cosso"
     return(cossoqrobj)
    }


ACV.lambda <- function(Gramat,time,status,mscale,basis.id,RS,cand.lambda0)
   {
   if(missing(cand.lambda0)) cand.lambda0=2^seq(-9,-18,-0.75)
   cand.lambda0=sort(cand.lambda0,decreasing=TRUE)
   n=length(time)
   p=length(mscale)

   tempBasisID=basis.id
   if(length(basis.id)>30)   tempBasisID=sort(sample(basis.id,30))

   Rtheta1=wsGram(Gramat[,tempBasisID,],mscale)

   Hess.FullNumer.unScale=array(NA,dim=c(length(tempBasisID),length(tempBasisID),n))
   for(i in 1:n)  Hess.FullNumer.unScale[,,i] =Rtheta1[i,]%*%t(Rtheta1[i,])

   ACV=matrix(NA,nrow=length(cand.lambda0),ncol=2)
   for(j in 1:length(cand.lambda0))
      {
      tempCox=sspline.Cox(Gramat[,tempBasisID,],Gramat[tempBasisID,tempBasisID,],time,status,mscale,tempBasisID,cand.lambda0[j],RS,Hess.FullNumer.unScale)
      ACV[j,1]=PartialLik(time,status,RS,tempCox$fit)
      ACV[j,2]=ACV[j,1]+sum(status==1)/n^2*( sum(diag(tempCox$UHU))/(n-1) - sum(tempCox$UHU)/(n^2-n) )
      }
   acv=ACV[,2]
   locMinid <- which((acv[-c(length(acv),length(acv)-1)]> acv[-c(1,length(acv))])*(acv[-c(1,length(acv))]<acv[-c(1:2)])==TRUE)+1
   locMaxid <- which((acv[-c(length(acv),length(acv)-1)]< acv[-c(1,length(acv))])*(acv[-c(1,length(acv))]>acv[-c(1:2)])==TRUE)+1
   locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(acv))]
   opt.lambda0=cand.lambda0[which.min(acv)]
   if(length(locMinid)>0)   opt.lambda0=cand.lambda0[ locMinid[which.min(acv[locMinid[1:length(locMinid)]])] ]

   return(opt.lambda0)
   }


tune.cosso.cox <- function(object,plot.it=TRUE,parallel=FALSE,cpus=1)
   {
    n <- nrow(object$x)
    nbasis <- length(object$basis.id)
    d <- length(object$wt)

    origMgrid=object$tune$Mgrid
    refinePt=which(apply(object$tune$L2norm<1e-8,1,sum)[-length(origMgrid)]-apply(object$tune$L2norm<1e-8,1,sum)[-1]>1)
    if(length(refinePt)==0)
       {   extMgrid <- c(apply(cbind(origMgrid[-1],origMgrid[-length(origMgrid)]),1,mean)[2:ceiling(length(origMgrid)/3+1)],max(origMgrid)+0.5) }
    else
       {   extMgrid <- c(as.numeric(apply(cbind(origMgrid[refinePt],origMgrid[refinePt+1]),1,quantile,c(.3,.6))),max(origMgrid)+0.5)  }

    if(length(extMgrid)<=5)  parallel=FALSE
    
    extACV <- rep(NA,length(extMgrid))
    extL2normMat <- matrix(NA,ncol=d,nrow=length(extMgrid))

     sfInit(parallel=parallel,cpus=cpus)
     if(parallel)
        {
        sfLibrary("quadprog",character.only=TRUE);sfLibrary("glmnet",character.only=TRUE)
        sfLibrary("cosso",character.only=TRUE)
        }
    tempcoefs <- sfClusterApplyLB(extMgrid,twostep.Cox,Gramat1=object$Kmat[,object$basis.id,],Gramat2=object$Kmat[object$basis.id,object$basis.id,],time=object$time,status=object$status,wt=object$wt,basis.id=object$basis.id,RS=object$RiskSet,lambda0=object$tune$OptLam)
    sfStop()
    for(jj in 1:length(extMgrid))
        {
        tempObj=tempcoefs[[jj]]
        for(j in 1:d)  extL2normMat[jj,j]=sqrt(mean((tempObj$theta[j]/object$wt[j]^2*object$Kmat[,object$basis.id,j]%*%tempObj$coef)^2))
        extACV[jj] <- PartialLik(object$time,object$status,object$RiskSet,tempObj$fit)+sum(object$status==1)/n^2*( sum(diag(tempObj$UHU))/(n-1) - sum(tempObj$UHU)/(n^2-n) )
        }

     Mgrid=c(origMgrid,extMgrid)
     ACV=c(object$tune$ACV,extACV)[order(Mgrid)]
     L2norm=rbind(object$tune$L2norm,extL2normMat)[order(Mgrid),]
     Mgrid=Mgrid[order( Mgrid)]

     locMinid <- which((ACV[-c(length(ACV),length(ACV)-1)]> ACV[-c(1,length(ACV))])*(ACV[-c(1,length(ACV))]<ACV[-c(1:2)])==TRUE)+1
     locMaxid <- which((ACV[-c(length(ACV),length(ACV)-1)]< ACV[-c(1,length(ACV))])*(ACV[-c(1,length(ACV))]>ACV[-c(1:2)])==TRUE)+1
     locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(ACV))]
     bestM=Mgrid[which.min(ACV)]
     if(length(locMinid)>0)   bestM=Mgrid[ locMinid[which.min(ACV[locMinid[1:length(locMinid)]])] ]
     tuneobj=list(OptM=bestM,OptLam=object$tune$OptLam,Mgrid=Mgrid,ACV=ACV,L2norm=L2norm)


     if(plot.it)
       {
       par(mfcol=c(1,2))
       plotid=complete.cases(cbind(tuneobj$Mgrid,tuneobj$ACV))
       plot(tuneobj$Mgrid[plotid],tuneobj$ACV[plotid],type="l",lwd=1.5,xlab="M",ylab="ACV")
       abline(v=tuneobj$OptM,lty=2,col=2);axis(3,tuneobj$OptM)
       matplot(tuneobj$Mgrid,tuneobj$L2norm,type="l",lty=1,col=c(1,rainbow(d-1)),xlab="M",ylab=expression(L[2]-norm))
       abline(v=tuneobj$OptM,lty=2,col=2);axis(3,tuneobj$OptM)
       axis(4,at=tuneobj$L2norm[length(tuneobj$Mgrid),],labels=1:d,cex=.3,las=2)
       }
     return(tuneobj)
    }

#---- Solve SS-ANOVA Cox and then a Quadratic Programming ----#
twostep.Cox <- function(Gramat1,Gramat2,time,status,wt,basis.id,RS,lambda0,M)
  {
  n=length(time)
  p=length(wt)
  #---- Step 1.2 ----#
  ssCox =sspline.Cox(Gramat1,Gramat2,time,status,    rep(1,p)/wt^2,basis.id,lambda0,RS)
  init.Theta=ssCox$L2norm*M/sum(ssCox$L2norm)
  #---- Step 2.1 ----#
   garCox=garrote.Cox(Gramat1,Gramat2,time,status,wt,basis.id,lambda0,M,ssCox$coef,init.Theta,RS)
  #---- Step 2.2 ----#
   ssCox =sspline.Cox(Gramat1,Gramat2,time,status,garCox$theta/wt^2,basis.id,lambda0,RS)
  obj=c(ssCox,list(theta=garCox$theta))

  return(obj)
  }

#---- Solve a SS-ANOVA Cox Problem ----#
sspline.Cox <- function(Gramat1,Gramat2,time,status,mscale,basis.id,lambda0,RS,Hess.FullNumer.unScale)
  {
  n=length(time)
  p=length(mscale)
  Rtheta1=wsGram(Gramat1,mscale)
  Rtheta2=wsGram(Gramat2,mscale)

  EigRtheta2=eigen(Rtheta2)
  if(min(EigRtheta2$value)<0)
     {
     Rtheta2=Rtheta2+max(1e-7,1.5*abs(min(EigRtheta2$value)))*diag(length(basis.id))
     EigRtheta2=eigen(Rtheta2)
     }
  pseudoX=Rtheta1%*%EigRtheta2$vectors%*%diag(sqrt(1/EigRtheta2$values))
  ssCox.en=glmnet(pseudoX,cbind(time=time,status=status),family="cox",lambda=c(lambda0/2,lambda0),alpha=0,standardize = FALSE)
  init.C=as.numeric( EigRtheta2$vectors%*%diag(sqrt(1/EigRtheta2$values))%*%ssCox.en$beta[,1] )

  #---- One-Step Update ----#
  f.old=Rtheta1%*%init.C
  GH=gradient.Hessian.C(init.C,Gramat1,Gramat2,time,status,mscale,lambda0,RS,Hess.FullNumer.unScale)
  new.C=My_solve(GH$H,GH$H%*%init.C-GH$G)

  #-------------------------#
  L2norm=rep(NA,p)
  for(j in 1:p)   L2norm[j]=sqrt(mean((mscale[j]*Gramat1[,,j]%*%new.C)^2))
  fit=Rtheta1%*%new.C

  UHU=Rtheta1%*%My_solve(GH$H,t(Rtheta1))
  ssCoxObj=list(coefs=new.C,fit=fit,L2norm=L2norm,UHU=UHU)
  return(ssCoxObj)
  }

gradient.Hessian.C=function(initC,Gramat1,Gramat2,time,status,mscale,lambda0,riskset,Hess.FullNumer.unScale)
  {
  n=length(time)
  tie.size=as.numeric( table(time[status==1]) ) 

  Rtheta1=wsGram(Gramat1,mscale)
  Rtheta2=wsGram(Gramat2,mscale)
  if(min(eigen(Rtheta2)$value)<0)  Rtheta2=Rtheta2+1e-8*diag(nrow(Rtheta2))
  eta=Rtheta1%*%initC

  if(missing(Hess.FullNumer.unScale))
     {
     Hess.FullNumer.unScale=array(NA,dim=c(length(initC),length(initC),n))
     for(i in 1:n)  Hess.FullNumer.unScale[,,i] =Rtheta1[i,]%*%t(Rtheta1[i,])
     }

  Grad.Term1=-t(Rtheta1)%*%status/n
  Grad.Term2=matrix(NA,ncol=length(riskset),nrow=length(initC))
  Grad.Term3=2*lambda0*Rtheta2%*%initC

  Grad.FullNumer=t(Rtheta1)%*%diag(as.numeric(exp(eta)))   
  Grad.FullDenom=Hess.FullDenom=exp(eta)                   

  Hess.FullNumer =Hess.FullNumer.unScale*array( rep( exp(eta),each=length(initC)^2 ), dim=c(length(initC),length(initC),n))
  Hess.Term1=Hess.Term2=array(NA,dim=c(length(initC),length(initC),length(riskset)))

  k=1
  tempSum.exp.eta=sum( exp(eta[ riskset[[k]] ]) )
  temp.Gradient.numer=apply(Grad.FullNumer[, riskset[[k]] ],1     ,sum)
  temp.Hessian.numer =apply(Hess.FullNumer[,,riskset[[k]] ],c(1,2),sum)

  Grad.Term2[,k] =tie.size[k]*temp.Gradient.numer/tempSum.exp.eta
  Hess.Term1[,,k]=temp.Hessian.numer /tempSum.exp.eta
  Hess.Term2[,,k]=1/tie.size[k]*Grad.Term2[,k]%*%t(Grad.Term2[,k])

  for(k in 2:length(riskset))
     {
     excludeID=riskset[[k-1]][ !riskset[[k-1]]%in%riskset[[k]] ]

     tempSum.exp.eta=tempSum.exp.eta-sum(exp(eta[excludeID]))
     if(length(excludeID)>1)
        {
        temp.Gradient.numer=temp.Gradient.numer-apply(Grad.FullNumer[, excludeID],1     ,sum)
        temp.Hessian.numer =temp.Hessian.numer -apply(Hess.FullNumer[,,excludeID],c(1,2),sum)
        }
     else
        {
        temp.Gradient.numer=temp.Gradient.numer-      Grad.FullNumer[, excludeID]
        temp.Hessian.numer =temp.Hessian.numer -      Hess.FullNumer[,,excludeID]
        }

      Grad.Term2[,k] =tie.size[k]*temp.Gradient.numer/tempSum.exp.eta
      Hess.Term1[,,k]=temp.Hessian.numer /tempSum.exp.eta
      Hess.Term2[,,k]=1/tie.size[k]*Grad.Term2[,k]%*%t(Grad.Term2[,k])
     }
  Grad.Term2=apply(Grad.Term2,1,sum)/n

  Gradient=Grad.Term1+Grad.Term2+Grad.Term3
  Hessian =apply(Hess.Term1,c(1,2),sum)/n-apply(Hess.Term2,c(1,2),sum)/n+2*lambda0*Rtheta2

  return(list(Gradient=Gradient,Hessian=Hessian))
  }


#---- Solve a Quadratic Programming Problem ------#
garrote.Cox <- function(Gramat1,Gramat2,time,status,wt,basis.id,lambda0,M,init.C,init.Theta,RS)
  {
  n=length(time)
  p=length(wt)

  if(missing(init.Theta))
     {
     L2norm=rep(NA,p)
     for(j in 1:p)   L2norm[j]=sqrt(mean((1/wt[j]^2*Gramat1[,,j]%*%init.C)^2))
     init.Theta=L2norm*M/sum(L2norm)
     }

  G1=matrix(NA,ncol=p,nrow=n)
  for(j in 1:p)     G1[,j]=1/wt[j]^2*Gramat1[,,j]%*%init.C
  G2=G1[basis.id,]

  Hess.FullNumer.unScale=array(NA,dim=c(length(init.Theta),length(init.Theta),n))
  for(i in 1:n)  Hess.FullNumer.unScale[,,i] =G1[i,]%*%t(G1[i,])

  loop=0
  iter.diff=Inf
  old.Theta=init.Theta
  while(loop<15 & iter.diff>1e-4)
      {
      loop=loop+1
      GH=gradient.Hessian.Theta(old.Theta,init.C,G1,G2,lambda0,M,time,status,RS,Hess.FullNumer.unScale)
      if(min(eigen(GH$H)$value)<0)   GH$H=GH$H+max( 1e-7,1.5*abs(min(eigen(GH$H)$value)) )*diag(length(init.Theta))
      dvec=-(GH$G-GH$H%*%old.Theta)
      Amat=t(rbind(diag(p),rep(-1,p)))
      bvec=c(rep(0,p),-M)
      new.Theta=My_solve.QP(GH$H,dvec,Amat,bvec)
      new.Theta[new.Theta<1e-7]=0
      iter.diff=mean(abs(new.Theta-old.Theta))
      old.Theta=new.Theta
      }
  return(list(coefs=init.C,theta=new.Theta))
  }

gradient.Hessian.Theta=function(initTheta,initC,G1,G2,lambda0,M,time,status,riskset,Hess.FullNumer.unScale)
  {
  n=length(time)
  p=length(initTheta)
  tie.size=as.numeric( table(time[status==1]) ) 
  eta=G1%*%initTheta

  Grad.Term1=-t(G1)%*%status/n
  Grad.Term2=matrix(NA,ncol=length(riskset),nrow=p)
  Grad.Term3=lambda0*t(G2)%*%initC

  Grad.FullNumer=t(G1)%*%diag(as.numeric(exp(eta)))   
  Grad.FullDenom=Hess.FullDenom=exp(eta)              

  Hess.FullNumer =Hess.FullNumer.unScale*array( rep( exp(eta),each=p^2 ), dim=c(p,p,n))
  Hess.Term1=Hess.Term2=array(NA,dim=c(p,p,length(riskset)))

  k=1
     tempSum.exp.eta=sum( exp( eta[ riskset[[k]] ] ) )
     tempGradient.numer=apply( Grad.FullNumer[, riskset[[k]] ],1     ,sum)
     tempHessian.numer =apply(Hess.FullNumer[,, riskset[[k]] ],c(1,2),sum)

     Grad.Term2[,k] =tie.size[k]*tempGradient.numer/tempSum.exp.eta
     Hess.Term1[,,k]=            tempHessian.numer /tempSum.exp.eta
     Hess.Term2[,,k]=1/tie.size[k]*Grad.Term2[,k]%*%t(Grad.Term2[,k])

  for(k in 2:length(riskset))
     {
     excludeID=riskset[[k-1]][! riskset[[k-1]]%in%riskset[[k]] ]
     tempSum.exp.eta=tempSum.exp.eta-sum( exp(eta[excludeID]) )

     if(length(excludeID)>1)
        {
        tempGradient.numer=tempGradient.numer-apply(Grad.FullNumer[, excludeID],1     ,sum)
        tempHessian.numer =tempHessian.numer -apply(Hess.FullNumer[,,excludeID],c(1,2),sum)
        }
     else
        {
        tempGradient.numer=tempGradient.numer-Grad.FullNumer[, excludeID]
        tempHessian.numer =tempHessian.numer -Hess.FullNumer[,,excludeID]
        }
     Grad.Term2[,k] =tie.size[k]*tempGradient.numer/tempSum.exp.eta
     Hess.Term1[,,k]=            tempHessian.numer /tempSum.exp.eta
     Hess.Term2[,,k]=1/tie.size[k]*Grad.Term2[,k]%*%t(Grad.Term2[,k])
     }
  Grad.Term2=apply(Grad.Term2,1,sum)/n

  Gradient=Grad.Term1+Grad.Term2+Grad.Term3
  Hessian =apply(Hess.Term1,c(1,2),sum)/n-apply(Hess.Term2,c(1,2),sum)/n

  return(list(Gradient=Gradient,Hessian=Hessian))
  }

SSANOVAwt.cox <- function(x,time,status,mscale=rep(1,ncol(x)))
  { 
    n <- length(time)
    d <- ncol(x)
    Gramat <- bigGram(x,x)
    basis.id <- sort(sample(which(status==1),min(30+ceiling(n/100),sum(status==1)-5)))
    RS <- RiskSet(time,status)
    optLambda <- ACV.lambda(Gramat,time,status,mscale,basis.id,RS,2^seq(-8,-20,-1))
    coxObj <- sspline.Cox(Gramat[,basis.id,],Gramat[basis.id,basis.id,],time,status,mscale,basis.id,optLambda,RS)
    L2norm <- rep(NA,d)
    for (j in 1:d)      L2norm[j] <- sqrt(mean((mscale[j]*Gramat[,basis.id,j]%*%coxObj$coef)^2))
  return(1/L2norm)
  }


PartialLik=function(time,status,RS,fhat)
   {
   pl=rep(NA,length(RS))
   eventtime=unique(time[status==1])
   tie.size=as.numeric(table(time[status==1]))
   for(k in 1:length(RS))
      {
      failid=which(time==eventtime[k])
      pl[k]=sum(fhat[failid])-tie.size[k]*log( sum(exp( fhat[RS[[k]]]) ) )
      }
   return(-sum(pl)/length(time))
   }



RiskSet <-  function(time,status)
  {
  eventTime=sort(unique(time[status==1]))
  RiskSet=list()
  for(k in 1:length(eventTime))
     {
     RiskSet=c(RiskSet,list(which(time>=eventTime[k])))
     }
  return(RiskSet)
  }

