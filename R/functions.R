
#' Make a model of vp ~ confining + deviatoric stress for mapping
#' @param df
#' @return a linear model of vp
#' @export
fit_vp_model <- function(df){
  lm(vp ~ confining_pressure + deviatoric_stress, data = df)
}

#' Make a model of vs ~ confining + deviatoric stress for mapping
#' @param df
#' @return a linear model of vs
#' @export
fit_vs_model <- function(df){
  lm(vs ~ confining_pressure + deviatoric_stress, data = df)
}

#' Root Mean Square Error of a model
#' @param model
#' @return scalar RMSE
#' @export
rmse <- function(model){
  sqrt(sum(model$residuals^2)/length(model$residuals))
}




