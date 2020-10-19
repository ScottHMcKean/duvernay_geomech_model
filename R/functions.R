
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

#' Plot partial dependence plot
#'
#' @param features feature vector (char vec)
#' @param target name of target (char)
#' @param predictor IML predictor
#' @param prefix named prefix of output directory
#' @return TRUE
#' @export
plot_pdp <- function(features, target, predictor, prefix){
  plist <- map(.x = 1:length(features), .f = make_pdp_plot,
               predictor = predictor, feats = features)

  ggsave(file = paste0('../output/',prefix,"/pdpplot.pdf"),
         arrangeGrob(grobs = plist, ncol = 3, left = target),
         width = 12, height = 6, dpi=600)

  TRUE
}

#' Plot importance and interation plot
#'
#' @param predictor IML predictor
#' @param prefix named prefix of output directory
#' @param predictor number of repetitions for permutation feature importance
#' @return TRUE
#' @export
plot_importance_interaction <- function(predictor, prefix, repetitions=500){
  importance = FeatureImp$new(predictor, loss='mae', n.repetitions = repetitions)

  p1 = ggplot(importance$results) +
    geom_col(aes(x = feature, y = importance)) +
    labs(x = 'Feature', y = 'Permutation Importance') +
    coord_flip() +
    theme_minimal()

  interact = Interaction$new(predictor)

  p2 = ggplot(as.data.frame(interact$results)) +
    geom_col(aes(x = .feature, y = .interaction)) +
    labs(x = '', y = 'Interaction') +
    coord_flip() +
    theme_minimal()

  ggsave(file = paste0('../output/',prefix,"/imp_int_plot.pdf"),
         arrangeGrob(grobs = list(p1,p2), ncol=1),
         width = 6, height = 6, dpi=600)

  TRUE
}

#' Plot residuals for train and test set
#'
#' @param train_predict output from predict(model, newdata=train)
#' @param test_predict output from predict(model, newdata=test)
#' @param prefix named prefix of output directory
#' @return TRUE
#' @export
plot_residuals <- function(train_predict, test_predict, prefix){
  p1 = plotResiduals(train_predict,loess.smooth=FALSE) +
    geom_abline(slope=1,linetype = "dashed") +
    ggtitle("Training Set") +
    theme_minimal()

  p2 = plotResiduals(test_predict,loess.smooth=FALSE) +
    geom_abline(slope=1,linetype = "dashed") +
    ggtitle("Test Set") +
    theme_minimal()

  ggsave(paste0('../output/',prefix,"/residual_plot.pdf"),
         arrangeGrob(grobs = list(p1,p2), nrow=1),
         width = 6, height = 6)
}

#' Wrapper Function for model generation
#'
#' @param target target for model (character)
#' @param features features for model (character vector)
#' @param model_type MLR embedded model type
#' @param pset MLR parameter set for model type
#' @param prefix named prefix of output directory
#' @param seed random seem for resampling
#' @param train_ratio ratio of train data to all data
#' @param mbo_budget budget for bayesian optimization
#' @param resample_iters number of iterations for resampling performance analysis
#' @param importance_reps number of iterations for permutation importance
#' @return TRUE
#' @export
generate_model <- function(
  target, features, model_type, pset, prefix, seed=1, train_ratio=0.9,
  mbo_budget= 75L, resample_iters=1000, importance_reps=50){
  set.seed(seed)
  ml_df = mldf_nores %>% select(features, target)

  train_rows = sample(nrow(ml_df)*train_ratio)
  train = ml_df[train_rows,]
  test = ml_df[-train_rows,]

  learner = makeLearner(model_type)

  train_task = makeRegrTask(data = train, target = target)
  test_task = makeRegrTask(data = test, target = target)
  all_task = makeRegrTask(data = ml_df, target = target)

  meas = list(mlr::mae, mlr::rmse)

  tune_res = tuneParams(
    learner, train_task, resampling=cv5, par.set=pset,
    control=makeTuneControlMBO(budget = mbo_budget),
    measures=meas
  )

  dir.create(paste0('../output/',prefix), showWarnings = FALSE)
  write_rds(tune_res, paste0('../output/',prefix,"/tune_results.rds"))

  # set hyperparameters
  tuned_learner = setHyperPars(learner = learner, par.vals = tune_res$x)

  # train final model for performance and interpretation
  model = train(tuned_learner, train_task)
  test_predict = predict(model, newdata=test)
  train_predict = predict(model, newdata=train)
  predictor = Predictor$new(model, data = ml_df)

  # monte-carlo performance measures
  resample = mlr::resample(
    tuned_learner, all_task,
    makeResampleDesc("Subsample", iters=resample_iters, split=4/5, predict='both'),
    measures = meas,
    show.info = FALSE
  )

  model_res = get_resample_regr_res(resample)
  model_res$train_perf = performance(predict(model, newdata=train), measures = meas)
  model_res$test_perf = performance(predict(model, newdata=test), measures = meas)
  model_res$model = model
  model_res$tuned_learner = tuned_learner

  sink(paste0('../output/',prefix,"/model_results.txt"))
  print(model_res)
  sink()

  write_rds(model_res, paste0('../output/',prefix,"/resample_results.rds"))

  # PLOT RESIDUALS
  plot_residuals(train_predict, test_predict, prefix)

  # IMPORTANCE & INTERACTION PLOTS
  plot_importance_interaction(predictor, prefix, repetitions=importance_reps)

  # PARTIAL DEPENDENCE PLOT
  plot_pdp(features, target, predictor, prefix)
}

