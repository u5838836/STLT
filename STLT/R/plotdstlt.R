#' @title Plot the Dyanmic Smooth Threshold Life Table
#'
#' @description Plot fitted DSTLT lines for the cohorts used to fit the model along with corresponding observed qx
#'
#' @param x a model object of class dstlt
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
#' dstlt.mod<-dstlt(cbind(60:100,60:100,60:100),cbind(seq(0.1,0.5,0.01),seq(0.1,0.5,0.01),seq(0.1,0.5,0.01)),endN=95)
#' plot(dstlt.mod)
#'
#' @export

plot.dstlt=function(x,main='',sub='',...){
  a=x$coefficients$a
  b=x$coefficients$b
  thet=x$coefficients$theta
  gam=x$coefficients$gamma
  N=x$coefficients$N
  start=x$Start
  taus=x$Taus
  omega=N-thet/gam
  qxs=x$qxs
  periods=x$periods

  if (abs(gam)<0.001) {
    # gam = 0
    X=seq(start,120,0.01)
    B=exp(a+b*0)
    C=1/(thet*B)^(1/N)
    gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
    gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
    paretoqx=(exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)-exp(-(X[((N-(start-1))*100+2):(5500+1)]-N)/thet))/exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)
    plot(X[1:(5400+1)],c(gompqx1,gompqx2,paretoqx),type="n",xlab="Age", ylab="qx",main=main, col='green',ylim=c(0,1),sub=sub)

    for (t in 1:periods) {
      B=exp(a+b*t)
      C=1/(thet*B)^(1/N)
      gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
      gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
      paretoqx=(exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)-exp(-(X[((N-(start-1))*100+2):(5500+1)]-N)/thet))/exp(-(X[((N-start)*100+2):(5400+1)]-N)/thet)
      lines(X[1:(5400+1)],c(gompqx,paretoqx),type="l",xlab="Age", ylab="qx",main=main, col=rainbow(periods)[t],xlim=c(60,120),ylim=c(0,1))
      points(start:taus[t],qxs[1:(which(is.na(qxs[,t]))[1]-1),t],col=rainbow(periods)[t])
    }
    #legend(start,0.7,#legend = 1:periods, col=rainbow(periods),pch=15)
  } else {
    if (gam<0) {
      #if gam<0
      X=seq(start,omega+1,0.01)
      B=exp(a+b*0)
      C=1/(thet*B)^(1/N)
      gompqx=(exp(-B/log(C)*(C^X[1:((N-start)*100+1)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1)]-1))
      paretoqx=((1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(100*((floor(100*(N-thet/gam))/100)-(start-1))+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)
      plot(X[1:(100*((floor(100*(N-thet/gam))/100)-start)+1)],c(gompqx,paretoqx),type="n",xlab="Age", ylab="qx",main=main,sub=sub)

      for (t in 1:periods) {
        B=exp(a+b*t)
        C=1/(thet*B)^(1/N)
        gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
        gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
        paretoqx=((1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(100*((floor(100*(N-thet/gam))/100)-(start-1))+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(100*((floor(100*(N-thet/gam))/100)-start)+1)]-N)/thet)^(-1/gam)
        lines(X[1:(100*((floor(100*(N-thet/gam))/100)-start)+1)],c(gompqx1,gompqx2,paretoqx),type="l",col=rainbow(periods)[t])
        points(start:taus[t],qxs[1:(which(is.na(qxs[,t]))[1]-1),t],col=rainbow(periods)[t])
      }
      #legend(start,0.7,#legend = 1:periods, col=rainbow(periods),pch=15)
    } else {
      # gam > 0

      X=seq(start,120,0.01)
      B=exp(a+b*0)
      C=1/(thet*B)^(1/N)
      gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
      gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
      paretoqx=((1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(5500+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)
      plot(X[1:(5400+1)],c(gompqx,paretoqx),type="n",col='blue',main=main,xlab="Age", ylab="qx",ylim=c(0,1),sub=sub)

      for (t in 1:periods) {
        B=exp(a+b*t)
        C=1/(thet*B)^(1/N)
        gompqx1=(exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))-exp(-B/log(C)*(C^X[101:((N-(start-1))*100+1-100)]-1)))/exp(-B/log(C)*(C^X[1:((N-start)*100+1-100)]-1))
        gompqx2=(-1+exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*(X[((N-start)*100+2):((N-start)*100+1+100)]-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^X[((N-start)*100+2-100):((N-start)*100+1)]-1))
        paretoqx=((1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)-(1+gam*(X[((N-(start-1))*100+2):(5500+1)]-N)/thet)^(-1/gam))/(1+gam*(X[((N-start)*100+2):(5400+1)]-N)/thet)^(-1/gam)
        lines(X[1:(5400+1)],c(gompqx1,gompqx2,paretoqx),type="l", col=rainbow(periods)[t], main=main,xlab="Age", ylab="qx")
        points(x,qx[1:(which(is.na(qx))[1]-1)])
      }
      #legend(start,0.7,#legend = 1:periods, col=rainbow(periods),pch=15)
    }
  }
}
