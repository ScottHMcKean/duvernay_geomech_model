
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

#' Calculate the bias and variance for a resample object on a regression
#'
#' @param resample resample object from mlr::resample, predict must be 'both'
#' @return list of train/test bias and variance + resample object for rds stash
#' @export
get_resample_regr_res = function(resample){
  train_pred = resample$pred$data %>%
    filter(set=='train') %>%
    group_by(id) %>%
    summarise(mean_pred = mean(response))

  train = resample$pred$data %>%
    filter(set=='train') %>%
    merge(train_pred, on='id') %>%
    mutate(ind_var = (response-mean_pred)^2,
           ind_bias = (response-truth)^2) %>%
    summarize(mean_bias = mean(ind_bias),
              mean_var = mean(ind_var)) %>%
    as.list()

  test_pred = resample$pred$data %>%
    filter(set=='test') %>%
    group_by(id) %>%
    summarise(mean_pred = mean(response))

  test = resample$pred$data %>%
    filter(set=='test') %>%
    merge(test_pred, on='id') %>%
    mutate(ind_var = (response-mean_pred)^2,
           ind_bias = (response-truth)^2) %>%
    summarize(mean_bias = mean(ind_bias),
              mean_var = mean(ind_var)) %>%
    as.list()

  return(
    list(train_mean_bias = train$mean_bias,
         train_mean_variance=train$mean_var,
         test_mean_bias=test$mean_bias,
         test_mean_variance=test$mean_var,
         resample_obj = resample)
  )
}

#' Make a partial dependence plot of all features
#'
#' @param i index for mapping
#' @param predictor IML predictor
#' @param feats feature vector (char vec)
#' @return individual pdp plot
#' @export
make_pdp_plot <- function(i, predictor, feats){
  pdp <- FeatureEffect$new(predictor, method = 'pdp+ice',
                           feature = feats[i])

  plot(pdp) + theme_minimal() + theme(axis.title.y=element_blank())
}
