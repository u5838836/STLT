#' @title Predict using a fitted Smooth Threshold Life Table
#'
#' @description For a given STLT fit, this function computes the qx associated with any particular age
#'
#' @param object an STLT model object
#'
#' @param newdata Vector of ages for which qx are to be predicted
#'
#' @param ... Additional arguments passed to methods
#'
#' @return Vector of predicted qx under the fitted STLT model
#'
#' @examples
#' stlt.mod<-stlt(60:100,seq(0.1,0.5,0.01))
#' predict(stlt.mod,65:90)
#'
#' @export

#prediction
predict.stlt=function(object,newdata,...){
  B=object$coefficients$B
  C=object$coefficients$C
  gam=object$coefficients$gamma
  N=object$coefficients$N
  start=object$Start
  tau=object$Tau
  thet=1/(C^N*B)
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





