#' @title Predict using the Dyanmic Smooth Threshold Life Table
#'
#' @description Returns the qx under a fitted DSTLT model, for a given vector of ages and scalar time period
#'
#' @param object a model object of class dstlt
#'
#' @param newdata vector of ages at which qx are to be predicted
#'
#' @param t the time period for which qx are to be predicted.t=1 corresponds to the earliest year to which the DSTLT model was fit; t=2 to the second earliest, etc.
#'
#' @param ... Additional arguments passed to methods
#'
#' @return A vector of predicted qx for a given numeric time period under the fitted DSTLT model
#'
#' @examples
#' dstlt.mod<-dstlt(cbind(60:100,60:100,60:100),cbind(seq(0.1,0.5,0.01),seq(0.1,0.5,0.01),seq(0.1,0.5,0.01)),endN=95)
#' predict(dstlt.mod,65:90,2)
#'
#' @export

#prediction
predict.dstlt=function(object,newdata,t,...){
  a=object$coefficients$a
  b=object$coefficients$b
  thet=object$coefficients$theta
  gam=object$coefficients$gamma
  N=object$coefficients$N
  B=exp(a+b*t)
  C=1/(thet*B)^(1/N)
  start=object$Start
  taus=object$Taus
  omega=N-thet/gam
  qx=rep(0,1)

  for (i in 1:length(newdata)) {
    x=newdata[i]
    if (x<N) {
      qx[i]=(exp(-B/log(C)*(C^x-1))-exp(-B/log(C)*(C^(x+1)-1)))/exp(-B/log(C)*(C^x-1))
    } else {
      if (x<N+1) {
        qx[i]=(-1+exp(-B/log(C)*(C^x-1))+1-(exp(-B/log(C)*(C^N-1)))*(1+gam*((x+1)-N)/thet)^(-1/gam))/exp(-B/log(C)*(C^x-1))
      } else {
        if (x<(omega-1)) {
          qx[i]=((1+gam*(x-N)/thet)^(-1/gam)-(1+gam*((x+1)-N)/thet)^(-1/gam))/(1+gam*(x-N)/thet)^(-1/gam)
        } else {
          if (x<omega) {
            qx[i]=1
          } else {
            qx[i]=NA
          }
        }
      }
    }
  }
  return(qx)
}







