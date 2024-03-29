% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictdstlt.R
\name{predict.dstlt}
\alias{predict.dstlt}
\title{Predict using the Dyanmic Smooth Threshold Life Table}
\usage{
\method{predict}{dstlt}(object, newdata, t, ...)
}
\arguments{
\item{object}{a model object of class dstlt}

\item{newdata}{vector of ages at which qx are to be predicted}

\item{t}{the time period for which qx are to be predicted.t=1 corresponds to the earliest year to which the DSTLT model was fit; t=2 to the second earliest, etc.}
}
\value{
A vector of predicted qx for a given numeric time period under the fitted DSTLT model
}
\description{
Returns the qx under a fitted DSTLT model, for a given vector of ages and scalar time period
}
\examples{
dstlt.mod<-dstlt(cbind(60:100,60:100,60:100),cbind(seq(0.1,0.5,0.01),seq(0.1,0.5,0.01),seq(0.1,0.5,0.01)),endN=95)
predict(dstlt.mod,65:90,2)

}
