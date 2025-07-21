#' @title Plot the Smooth Threshold Life Table
#'
#' @description Plots the fitted threshold life table line and corresponding observed qx used for fitting the model.
#'
#' @param x An stlt model object
#'
#' @param main title of plot
#'
#' @param sub subtitle of plot
#'
#' @param ... Additional arguments passed to methods
#'
#' @return NULL
#'
#' @examples
#' stlt.mod<-stlt(60:100,seq(0.1,0.5,0.01))
#' plot(stlt.mod)
#'
#' @export


plot.stlt=function(x,main='',sub='',...){
  B=x$coefficients$B
  C=x$coefficients$C
  gam=x$coefficients$gamma
  N=x$coefficients$N
  start=x$Start
  tau=x$Tau
  thet=1/(C^N*B)
  omega=N-thet/gam
  qx=x$qx

  if (abs(gam)<0.001) {
    # gam = 0
    x=start:tau
    X=seq(start,120)
    gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
    gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
    paretoqx=(exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)-exp(-(X[((N-(start-1))*100+2):(5500+1)]-N)/thet))/exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)
    plot(X[1:(5400+1)],c(gompqx1,gompqx2,paretoqx),type="l",xlab="Age", ylab="qx", col='black',ylim=c(0,1),main=main,sub=sub)
    points(x,qx[1:(which(is.na(qx))[1]-1)])
  } else {
    if (gam<0) {
      #if gam<0
      x=start:tau
      X=seq(start,omega+1,0.01)
      gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
      gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
      paretoqx=((1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(100*((floor(100*(N-thet/gam))/100)-(start-1))+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)
      plot(X[1:(100*((floor(100*(N-thet/gam))/100)-start)+1)],c(gompqx1,gompqx2,paretoqx),type="l",col='black',xlab='Age',ylab='qx',main=main,sub=sub)
      points(x,qx[1:(which(is.na(qx))[1]-1)])

    } else {
      # gam > 0
      x=start:tau
      X=seq(start,120,0.01)
      gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
      gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
      paretoqx=((1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(5500+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)
      plot(X[1:(5400+1)],c(gompqx1,gompqx2,paretoqx),type="l",col='black',ylim=c(0,1),main=main,sub=sub)
      points(x,qx[1:(which(is.na(qx))[1]-1)])
    }
  }
}
