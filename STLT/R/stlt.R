#' @title Fit the Smooth Threshold Life Table
#'
#' @description This package fits the Smooth Threshold Life Table and Dyanmic Smooth Threshold Life Table, as in Huang et al., (2020). Fitted and predicted qx as well as their plots are provided. No right censoring is applied and all possible ages are considered when estimating the threshold age N.
#'
#' @param ages Vector of ages that are used to fit the STLT. There should not be missing data for any ages between the smallest and largest age.
#'
#' @param qx Vector of mortality rates that are used to fit the STLT, corresponding to the vector of ages. The qx vector need not be equal to one for any element; right censoring is accounted for in the likelihood function.
#'
#' @param startN The first potential threshold age N tested
#'
#' @param endN The last potential threshold age N tested
#'
#' @param censorAge The age at which to right-censor observations
#'
#' @param hessian Set to TRUE for standard errors
#'
#' @param radix The radix for the hypothetical cohort; affects standard errors
#'
#' @return A model object of class "stlt"
#'
#' @examples stlt(60:100,seq(0.1,0.5,0.01))
#'
#' @export stlt
#
# @export plot.stlt
#
# @export predict.stlt

stlt<-function(ages,qx,startN=0,endN=200,censorAge=NULL,hessian=FALSE,radix=50000)
{
  start=ages[1]
  if (!is.null(censorAge)) {
  qx=qx[1:which(ages==censorAge)]
  }
  qx=append(qx,NA)
  lstart=radix
  dx=rep(0,50)
  exposures=rep(0,50)
  for (i in 1:(which(is.na(qx))[1]-1)) {
    exposures[i]=radix
    dx[i]=qx[i]*radix
    radix=radix-dx[i]
  }
  ltau=radix
  tau=which(is.na(qx))[1]-2+start

  #####choosing N####

  lps=rep(0,1)
  for (N in max(start+1,startN):min(tau-1,endN)) {

    x=start:(N-1)
    x1=(start+1):(N)
    x2=N:tau
    x3=(N+1):(tau+1)

    log.lik <- function(theta){
      alpha <- theta[1]
      delta <- theta[2]
      gam <- theta[3]
      out <- sum(dx[1:(N-start)]*(-(exp(alpha))/(exp(delta))*((exp(exp(delta)))^x-1)+log(1-exp(-(exp(alpha))/(exp(delta))*(exp(exp(delta)))^x*((exp(exp(delta)))-1)))))+exposures[N-(start-1)]*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^N-1))-lstart*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^start-1))+sum(dx[(N-(start-1)):(which(is.na(qx))[1]-1)]*log((1+gam*((x2-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)-(1+gam*((x3-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)))+ltau*log((1+gam*((tau-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam))
      return(out)
    }

    theta.start <- c(-11,-3,2)
    out <- optim(theta.start, log.lik, hessian = FALSE, control = list(fnscale=-1), method="Nelder-Mead")
    beta.hat <- out$par
    beta.hat
    alpha=beta.hat[1]
    delta=beta.hat[2]
    B=exp(beta.hat[1])
    C=exp(exp(beta.hat[2]))
    gam=beta.hat[3]

    lps[N-max(start+1,startN)+1]=sum(dx[1:(N-start)]*(-(exp(alpha))/(exp(delta))*((exp(exp(delta)))^x-1)+log(1-exp(-(exp(alpha))/(exp(delta))*(exp(exp(delta)))^x*((exp(exp(delta)))-1)))))+exposures[N-(start-1)]*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^N-1))-lstart*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^start-1))+sum(dx[(N-(start-1)):(which(is.na(qx))[1]-1)]*log((1+gam*((x2-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)-(1+gam*((x3-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)))+ltau*log((1+gam*((tau-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam))
  }

  stlt=max(lps)
  which.max(lps)
  N=max(start+1,startN)+which.max(lps)-1
  NSE=(1/((lps[which.max(lps)]-lps[which.max(lps)-1])-(lps[which.max(lps)+1]-lps[which.max(lps)])))^(1/2)

  ####refit with optimal N####

  x=start:(N-1)
  x1=(start+1):(N)
  x2=N:tau
  x3=(N+1):(tau+1)

  log.lik <- function(theta){
    alpha <- theta[1]
    delta <- theta[2]
    gam <- theta[3]
    out <- sum(dx[1:(N-start)]*(-(exp(alpha))/(exp(delta))*((exp(exp(delta)))^x-1)+log(1-exp(-(exp(alpha))/(exp(delta))*(exp(exp(delta)))^x*((exp(exp(delta)))-1)))))+exposures[N-(start-1)]*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^N-1))-lstart*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^start-1))+sum(dx[(N-(start-1)):(which(is.na(qx))[1]-1)]*log((1+gam*((x2-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)-(1+gam*((x3-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)))+ltau*log((1+gam*((tau-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam))
    return(out)
  }
  if (hessian==TRUE) {
    theta.start <- c(-11,-3,2)
    out <- optim(theta.start, log.lik, hessian = TRUE, control = list(fnscale=-1), method="Nelder-Mead")
    beta.hat <- out$par
    beta.hat
    alpha=beta.hat[1]
    delta=beta.hat[2]
    B=exp(beta.hat[1])
    C=exp(exp(beta.hat[2]))
    gam=beta.hat[3]
    thet=1/(C^N*B)
    omega=N-thet/gam
    var.beta.hat <- diag(solve(-out$hessian))
    var.beta.hat
    se.beta.hat <- sqrt(var.beta.hat)
    se.beta.hat
    varcovar=solve(-out$hessian)
    omegavar=t(c(1/(gam*C^N*B),(N*log(C))/(gam*C^N*B),1/(gam^2*C^N*B)))%*%varcovar%*%c(1/(gam*C^N*B),(N*log(C))/(gam*C^N*B),1/(gam^2*C^N*B))
    omegase=sqrt(omegavar)

    returnlist=list(coefficients=list(B=B,C=C,gamma=gam,N=N),Omega=omega,Start=start,Tau=tau,SEs=list(B=exp(2*alpha)*se.beta.hat[1],C=exp(2*exp(delta))*exp(2*delta)*se.beta.hat[2],gam=se.beta.hat[3],Omega=omegase,N=NSE),qx=qx)
    class(returnlist)="stlt"
    return(returnlist)
  } else {
    theta.start <- c(-11,-3,2)
    out <- optim(theta.start, log.lik, hessian = FALSE, control = list(fnscale=-1), method="Nelder-Mead")
    beta.hat <- out$par
    beta.hat
    alpha=beta.hat[1]
    delta=beta.hat[2]
    B=exp(beta.hat[1])
    C=exp(exp(beta.hat[2]))
    gam=beta.hat[3]
    thet=1/(C^N*B)
    omega=N-thet/gam
    returnlist=list(coefficients=list(B=B,C=C,gamma=gam,N=N),Omega=omega,Start=start,Tau=tau,qx=qx)
    class(returnlist)="stlt"
    return(returnlist)
  }
}


# plot.stlt=function(x,...){
#   B=x$coefficients$B
#   C=x$coefficients$C
#   gam=x$coefficients$gamma
#   N=x$coefficients$N
#   start=x$Start
#   tau=x$Tau
#   thet=1/(C^N*B)
#   omega=N-thet/gam
#   qx=x$qx
#
#   if (abs(gam)<0.001) {
#     # gam = 0
#     x=start:tau
#     X=seq(start,120)
#     gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
#     gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
#     paretoqx=(exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)-exp(-(X[((N-(start-1))*100+2):(5500+1)]-N)/thet))/exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)
#     plot(X[1:(5400+1)],c(gompqx,paretoqx),type="l",xlab="Age", ylab="qx",main="STLT plot", col='blue',ylim=c(0,1))
#     points(x,qx[1:(which(is.na(qx))[1]-1)])
#   } else {
#     if (gam<0) {
#       #if gam<0
#       x=start:tau
#       X=seq(start,omega+1,0.01)
#       gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
#       gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
#       paretoqx=((1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(100*((floor(100*(N-thet/gam))/100)-(start-1))+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)
#       plot(X[1:(100*((floor(100*(N-thet/gam))/100)-start)+1)],c(gompqx1,gompqx2,paretoqx),type="l",col='blue',main='STLT plot',xlab='Age',ylab='qx')
#       points(x,qx[1:(which(is.na(qx))[1]-1)])
#
#     } else {
#       # gam > 0
#       x=start:tau
#       X=seq(start,120,0.01)
#       gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
#       gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
#       paretoqx=((1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(5500+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)
#       plot(X[1:(5400+1)],c(gompqx,paretoqx),type="l",col='red',main="STLT plot",ylim=c(0,1))
#       points(x,qx[1:(which(is.na(qx))[1]-1)])
#     }
#   }
# }
#
#
# #prediction
# predict.stlt=function(object,newdata,...){
#   B=object$coefficients$B
#   C=object$coefficients$C
#   gam=object$coefficients$gamma
#   N=object$coefficients$N
#   start=object$Start
#   tau=object$Tau
#   thet=1/(C^N*B)
#   omega=N-thet/gam
#   qx=rep(0,1)
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





