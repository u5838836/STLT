#' @title Get predictions from MortalityLaws models
#'
#' @description Extract predicted qx from models fit using the MortalityLaws package.
#'
#' @param x Vector of ages at the beginning of the age interval.
#'
#' @param qx Vector of observed mortality rates.
#'
#' @param law The law to get predictions from - see MortalityLaws package.
#'
#' @param pred_ages Ages at which to predict qx.
#'
#' @return Predicted qx at pred_ages from chosen law.
#'
#' @examples
#' qx <- seq(0.1,0.5,0.01)
#' get_qx(60:100, qx, "HP2", 80:105)
#'
#' @export get_qx

get_qx<-function(x, qx, law, pred_ages){
  fit = MortalityLaws::MortalityLaw(x = x, qx = qx, law = law)
  preds = predict(fit, x = pred_ages)
  if (law == 'HP2' | law == 'HP3' | law == 'HP4') {
    return(preds)
  } else {
    return(1-exp(-preds))
  }
}


