#' @title Plot fitted lines of any number of models
#'
#' @description Plot any number of models fit using the MortalityLaws package and the STLT.
#'
#' @param x Vector of ages at the beginning of the age interval.
#'
#' @param qx Vector of observed mortality rates.
#'
#' @param laws Vector of strings that indicate which laws to plot.
#'
#' @return Plot of fitted lines for chosen laws and STLT.
#'
#' @examples
#' plot_all(60:100,seq(0.1,0.5,0.01),c('HP2','gompertz'))
#'
#' @export plot_all

plot_all<-function(x, qx, laws){
  stlt_mod = stlt(ages = x, qx = qx)
  start=stlt_mod$Start
  plot_ages = seq(start,120,0.01)
  plot(stlt_mod)

  for (i in 1:length(laws)) {
    pred = get_qx(x=x,qx=qx,law=laws[i],pred_ages=plot_ages)
    lines(plot_ages,pred,col=i+1)
  }

  legend("topright", legend = c('STLT',laws), col = 1:(length(laws)+1), lty = 1)
}
