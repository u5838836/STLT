#' @title Fit the Dyanmic Smooth Threshold Life Table
#'
#' @description This package fits the Smooth Threshold Life Table and Dyanmic Smooth Threshold Life Table, as in Huang et al., (2020). Fitted and predicted qx as well as their plots are provided. No right censoring is applied and all possible ages are considered when estimating the threshold age N.
#'
#' @param ages Matrix of ages that are used to fit the DSTLT. Each column should represent one cohort. There should not be missing data for any ages between the smallest and largest age.
#'
#' @param qxs Matrix of mortality rates that are used to fit the DSTLT, corresponding to the matrix of ages. Each column should represent one cohort.The qx matrix need not be equal to one for any element; right censoring is accounted for in the likelihood function. Where the highest observed age at death is different across cohorts, elements in the qx matrix should be filled with NA for ages past the highest observed age for each cohort.
#'
#' @param startN The first potential threshold age N tested
#'
#' @param endN The last potential threshold age N tested
#'
#' @param censorAge The age at which to right-censor observations
#'
#' @param hessian set to TRUE for standard errors
#'
#' @return NULL
#'
#' @examples dstlt(cbind(60:100,60:100,60:100),cbind(seq(0.11,0.51,0.01),seq(0.105,0.505,0.01),seq(0.1,0.5,0.01)),endN=95)
#'
#' @export dstlt
#
# @export plot.dstlt
#
# @export predict.dstlt

dstlt<-function(ages,qxs,startN=80,endN=105,censorAge=NULL,hessian=FALSE)
{
  ages=rbind(ages,NA)
  if (!is.null(censorAge)) {
  qxs=qxs[1:which(ages[,1]==censorAge),]
  }
  qxs=rbind(qxs,NA)
  periods=ncol(ages)
  start=ages[1,1]
  ltaus=rep(0,periods)
  taus=rep(0,periods)
  dxs=matrix(nrow=200,ncol=periods)
  exposuress=matrix(nrow=200,ncol=periods)
  for (m in 1:periods) {
    taus[m]=which(is.na(qxs[,m]))[1]-2+start
    dx=rep(0,200)
    exposures=rep(0,200)
    radix=50000
    for (i in 1:(which(is.na(qxs[,m]))[1]-1)) {
      exposures[i]=radix
      dx[i]=qxs[i,m]*radix
      radix=radix-dx[i]
    }
    ltaus[m]=radix
    dxs[,m]=dx
    exposuress[,m]=exposures
  }

  #####choosing N####

  lps=rep(0,1)
  for (N in startN:endN) {

    x=start:(N-1)
    x1=x+1
    x2=matrix(nrow=200,ncol=periods)
    for (i in 1:periods) {
      x2[1:length(N:taus[i]),i]=N:taus[i]
    }
    x3=x2+1

    log.lik <- function(theta){
      a <- theta[1]
      b <- theta[2]
      thet <- theta[3]
      gam <- theta[4]
      timelike=rep(0,periods)
      for (t in 1:periods) {
        timelike[t] <- sum(dxs[1:(N-start),t]*(-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t))*(((thet*exp(a+b*t))^(-1/N))^x-1)+log(1-exp(-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t))*((thet*exp(a+b*t))^(-1/N))^x*(((thet*exp(a+b*t))^(-1/N))-1)))))+exposuress[N-(start-1),t]*((-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t)))*(((thet*exp(a+b*t))^(-1/N))^N-1))-exposuress[1,t]*((-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t)))*(((thet*exp(a+b*t))^(-1/N))^start-1))+sum(dxs[(N-(start-1)):(which(is.na(qxs[,t]))[1]-1),t]*log((1+gam*((x2[!is.na(x2[,t]),t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam)-(1+gam*((x3[!is.na(x3[,t]),t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam)))+ltaus[t]*log((1+gam*((taus[t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam))
      }
      out=sum(timelike)
      return(out)
    }

    theta.start <- c(-11,-2,3,2)
    out <- optim(theta.start, log.lik, hessian = FALSE, control = list(fnscale=-1), method="Nelder-Mead")
    beta.hat <- out$par
    beta.hat
    a=beta.hat[1]
    b=beta.hat[2]
    thet=beta.hat[3]
    gam=beta.hat[4]
    timelike=rep(0,periods)
    for (t in 1:periods) {
      timelike[t] <- sum(dxs[1:(N-start),t]*(-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t))*(((thet*exp(a+b*t))^(-1/N))^x-1)+log(1-exp(-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t))*((thet*exp(a+b*t))^(-1/N))^x*(((thet*exp(a+b*t))^(-1/N))-1)))))+exposuress[N-(start-1),t]*((-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t)))*(((thet*exp(a+b*t))^(-1/N))^N-1))-exposuress[1,t]*((-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t)))*(((thet*exp(a+b*t))^(-1/N))^start-1))+sum(dxs[(N-(start-1)):(which(is.na(qxs[,t]))[1]-1),t]*log((1+gam*((x2[!is.na(x2[,t]),t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam)-(1+gam*((x3[!is.na(x3[,t]),t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam)))+ltaus[t]*log((1+gam*((taus[t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam))
    }
    lps[N-startN+1]=sum(timelike)
  }

  DSTLTLP=max(lps)
  which.max(lps)
  N=startN+which.max(lps)-1
  NSE=(1/((lps[which.max(lps)]-lps[which.max(lps)-1])-(lps[which.max(lps)+1]-lps[which.max(lps)])))^(1/2)

  ####refit with optimal N####

  x=start:(N-1)
  x1=x+1
  x2=matrix(nrow=200,ncol=periods)
  for (i in 1:periods) {
    x2[1:length(N:taus[i]),i]=N:taus[i]
  }
  x3=x2+1

  log.lik <- function(theta){
    a <- theta[1]
    b <- theta[2]
    thet <- theta[3]
    gam <- theta[4]
    timelike=rep(0,periods)
    for (t in 1:periods) {
      timelike[t] <- sum(dxs[1:(N-start),t]*(-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t))*(((thet*exp(a+b*t))^(-1/N))^x-1)+log(1-exp(-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t))*((thet*exp(a+b*t))^(-1/N))^x*(((thet*exp(a+b*t))^(-1/N))-1)))))+exposuress[N-(start-1),t]*((-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t)))*(((thet*exp(a+b*t))^(-1/N))^N-1))-exposuress[1,t]*((-(exp(a+b*t))/((-1/N)*(log(thet)+a+b*t)))*(((thet*exp(a+b*t))^(-1/N))^start-1))+sum(dxs[(N-(start-1)):(which(is.na(qxs[,t]))[1]-1),t]*log((1+gam*((x2[!is.na(x2[,t]),t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam)-(1+gam*((x3[!is.na(x3[,t]),t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam)))+ltaus[t]*log((1+gam*((taus[t]-N)/(1/(((thet*exp(a+b*t))^(-1/N))^N*exp(a+b*t)))))^(-1/gam))
    }
    out=sum(timelike)
    return(out)
  }
  if (hessian==TRUE) {
    theta.start <- c(-11,-2,3,2)
    out <- optim(theta.start, log.lik, hessian = TRUE, control = list(fnscale=-1), method="Nelder-Mead")
    beta.hat <- out$par
    beta.hat
    a=beta.hat[1]
    b=beta.hat[2]
    thet=beta.hat[3]
    gam=beta.hat[4]
    omega=N-thet/gam
    var.beta.hat <- diag(solve(-out$hessian))
    se.beta.hat <- sqrt(var.beta.hat)
    returnlist=list(coefficients=list(a=a,b=b,theta=thet,gamma=gam,N=N),Omega=omega,Start=start,Taus=taus,SEs=list(a=se.beta.hat[1],b=se.beta.hat[2],theta=se.beta.hat[3],gamma=se.beta.hat[4],N=NSE),qxs=qxs,periods=periods)
    class(returnlist)="dstlt"
    return(returnlist)
  } else {
    theta.start <- c(-11,-2,3,2)
    out <- optim(theta.start, log.lik, hessian = FALSE, control = list(fnscale=-1), method="Nelder-Mead")
    beta.hat <- out$par
    beta.hat
    a=beta.hat[1]
    b=beta.hat[2]
    thet=beta.hat[3]
    gam=beta.hat[4]
    omega=N-thet/gam
    returnlist=list(coefficients=list(a=a,b=b,theta=thet,gamma=gam,N=N),Omega=omega,Start=start,Taus=taus,qxs=qxs,periods=periods)
    class(returnlist)="dstlt"
    return(returnlist)
  }
}


################ plot the DSTLT ####################

# plot.dstlt=function(x,...){
#   a=x$coefficients$a
#   b=x$coefficients$b
#   thet=x$coefficients$theta
#   gam=x$coefficients$gamma
#   N=x$coefficients$N
#   start=x$Start
#   taus=x$Taus
#   omega=N-thet/gam
#   qxs=x$qxs
#   periods=x$periods
#
#   if (abs(gam)<0.001) {
#     # gam = 0
#     X=seq(start,120,0.01)
#     B=exp(a+b*0)
#     C=1/(thet*B)^(1/N)
#     gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
#     gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
#     paretoqx=(exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)-exp(-(X[((N-(start-1))*100+2):(5500+1)]-N)/thet))/exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)
#     plot(X[1:(5400+1)],c(gompqx,paretoqx),type="n",xlab="Age", ylab="qx",main="DSTLT plot", col='blue',ylim=c(0,1))
#
#     for (t in 1:periods) {
#       B=exp(a+b*t)
#       C=1/(thet*B)^(1/N)
#       gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
#       gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
#       paretoqx=(exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)-exp(-(X[((N-(start-1))*100+2):(5500+1)]-N)/thet))/exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)
#       lines(X[1:(5400+1)],c(gompqx,paretoqx),type="l",xlab="Age", ylab="qx",main="DSTLT plot", col=rainbow(periods)[t],xlim=c(60,120),ylim=c(0,1))
#       points(start:taus[t],qxs[1:(which(is.na(qxs[,t]))[1]-1),t],col=rainbow(periods)[t])
#     }
#     legend(start,0.7,legend = 1:periods, col=rainbow(periods),pch=15)
#   } else {
#     if (gam<0) {
#       #if gam<0
#       X=seq(start,omega+1,0.01)
#       B=exp(a+b*0)
#       C=1/(thet*B)^(1/N)
#       gompqx=(exp(-B/log(C)*(C^X[1:((N-start)*100+1)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1)]-1))
#       paretoqx=((1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(100*((floor(100*(N-thet/gam))/100)-(start-1))+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)
#       plot(X[1:(100*((floor(100*(N-thet/gam))/100)-start)+1)],c(gompqx,paretoqx),type="n",xlab="Age", ylab="qx",main="DSTLT plot")
#
#       for (t in 1:periods) {
#         B=exp(a+b*t)
#         C=1/(thet*B)^(1/N)
#         gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
#         gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
#         paretoqx=((1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(100*((floor(100*(N-thet/gam))/100)-(start-1))+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)
#         lines(X[1:(100*((floor(100*(N-thet/gam))/100)-start)+1)],c(gompqx1,gompqx2,paretoqx),type="l",col=rainbow(periods)[t])
#         points(start:taus[t],qxs[1:(which(is.na(qxs[,t]))[1]-1),t],col=rainbow(periods)[t])
#       }
#       legend(start,0.7,legend = 1:periods, col=rainbow(periods),pch=15)
#     } else {
#       # gam > 0
#
#       X=seq(start,120,0.01)
#       B=exp(a+b*0)
#       C=1/(thet*B)^(1/N)
#       gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
#       gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
#       paretoqx=((1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(5500+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)
#       plot(X[1:(5400+1)],c(gompqx,paretoqx),type="n",col='blue',main="DSTLT plot",xlab="Age", ylab="qx",ylim=c(0,1))
#
#       for (t in 1:periods) {
#         B=exp(a+b*t)
#         C=1/(thet*B)^(1/N)
#         gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
#         gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
#         paretoqx=((1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(5500+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)
#         lines(X[1:(5400+1)],c(gompqx,paretoqx),type="l", col=rainbow(periods)[t], main="DSTLT plot",xlab="Age", ylab="qx")
#         points(x,qx[1:(which(is.na(qx))[1]-1)])
#       }
#       legend(start,0.7,legend = 1:periods, col=rainbow(periods),pch=15)
#     }
#   }
# }
#
#
#
# #prediction
# predict.dstlt=function(object,newdata,t,...){
#   a=object$coefficients$a
#   b=object$coefficients$b
#   thet=object$coefficients$theta
#   gam=object$coefficients$gamma
#   N=object$coefficients$N
#   B=exp(a+b*t)
#   C=1/(thet*B)^(1/N)
#   start=object$Start
#   taus=object$Taus
#   omega=N-thet/gam
#   qx=rep(0,1)
#
#   for (i in 1:length(newdata)) {
#     x=newdata[i]
#     if (x<N) {
#       qx[i]=(exp(-B/log(C)*(C^x-1))-exp(-B/log(C)*(C^(x+1)-1)))/exp(-B/log(C)*(C^x-1))
#     } else {
#       if (x<N+1) {
#         qx[i]=(-1+exp(-B/log(C)*(C^x-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*((x+1)-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^x-1))
#       } else {
#         if (x<(omega-1)) {
#           qx[i]=((1+gam*(x-N)/thet)^(-1/gam)-(1+gam*((x+1)-N)/thet)^(-1/gam))/(1+gam*(x-N)/thet)^(-1/gam)
#         } else {
#           if (x<omega) {
#             qx[i]=1
#           } else {
#             qx[i]=NA
#           }
#         }
#       }
#     }
#   }
#   return(qx)
# }







