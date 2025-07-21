#' @title Compare any number of models
#'
#' @description Compare any number of models fit using the MortalityLaws package and the STLT.
#'
#' @param x Vector of ages at the beginning of the age interval.
#'
#' @param qx Vector of observed mortality rates.
#'
#' @param laws Vector of strings that indicate which laws to compare.
#'
#' @param wx Vector of weights for the weighted prediction error measures. Usually the cohort size at each x. Can be a scalar indicating the cohort size at the first age - subsequent weights will then be calculated using the observed qx
#'
#' @param compare_ages Vector of ages to restrict the comparison to.
#'
#' @return Comparison of prediction accuracy measures for chosen laws and STLT.
#'
#' @examples
#' compare_all(60:100,seq(0.1,0.5,0.01),c('HP2','gompertz'),10000)
#'
#' @export compare_all

compare_all<-function(x,qx,laws,wx=NULL, compare_ages = NULL){
  if (is.null(compare_ages)) {
    compare_ages = x
  }

  stlt_mod = stlt(ages = x, qx = qx)
  stlt_preds = predict(stlt_mod, newdata = compare_ages)

  law_preds = matrix(nrow=length(compare_ages), ncol=length(laws))
  for (i in 1:length(laws)) {
    pred = get_qx(x=x,qx=qx,law=laws[i],pred_ages=compare_ages)
    law_preds[,i] = pred
  }
  colnames(law_preds) = laws

  all_preds = cbind(stlt_preds, law_preds)
  colnames(all_preds)[1] = "STLT"

  return(compare_methods(qx = qx, pred_qx = all_preds, wx = wx))
}


