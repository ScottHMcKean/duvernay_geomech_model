
#' Make a model of vp ~ confining + deviatoric stress for mapping
#' @param df machine learning dataframe
#' @return a linear model of vp
#' @export
fit_vp_model <- function(df){
  lm(vp ~ confining_pressure + deviatoric_stress, data = df)
}

#' Make a model of vs ~ confining + deviatoric stress for mapping
#' @param df machine learning dataframe
#' @return a linear model of vs
#' @export
fit_vs_model <- function(df){
  lm(vs ~ confining_pressure + deviatoric_stress, data = df)
}

#' Root Mean Square Error of a model
#' @param model model output with a residuals column / vector
#' @return scalar RMSE
#' @export
rmse <- function(model){
  sqrt(sum(model$residuals^2)/length(model$residuals))
}

#' Calculate the bias and variance for a resample object on a regression
#' @param resample resample object from mlr::resample, predict must be 'both'
#' @return list of train/test bias and variance + resample object for rds stash
#' @export
get_resample_regr_res = function(resample){
  train_pred = resample$pred$data %>%
    filter(set=='train') %>%
    group_by(id) %>%
    summarise(mean_pred = mean(response))

  train_res= resample$pred$data %>%
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

  test_res = resample$pred$data %>%
    filter(set=='test') %>%
    merge(test_pred, on='id') %>%
    mutate(ind_var = (response-mean_pred)^2,
           ind_bias = (response-truth)^2) %>%
    summarize(mean_bias = mean(ind_bias),
              mean_var = mean(ind_var)) %>%
    as.list()

  return(
    list(train_mean_bias = train_res$mean_bias,
         train_mean_variance= train_res$mean_var,
         test_mean_bias = test_res$mean_bias,
         test_mean_variance = test_res$mean_var,
         resample_obj = resample)
  )
}

#' Make a partial dependence plot of all features
#' @param i index for mapping
#' @param predictor IML predictor
#' @param feats feature vector (char vec)
#' @return individual pdp plot
#' @export
make_pdp_plot <- function(i, predictor, feats){
  pdp <- FeatureEffect$new(predictor, method = 'pdp+ice', center.at = 0,
                           feature = feats[i], grid.size=50)

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
  plist <- purrr::map(.x = 1:length(features), .f = make_pdp_plot,
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

  # write out importance results
  imp_int_res = importance$results %>%
    merge(., data.frame(interact$results), by.x='feature', by.y='.feature')

  write_csv(imp_int_res, paste0('../output/',prefix,"/imp_int_res.csv"))

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
generate_model <- function(data,
  target, features, model_type, pset, prefix, seed=1, train_ratio=0.9,
  mbo_budget= 75L, resample_iters=1000, importance_reps=50){
  set.seed(seed)

  ml_df = data %>%
    dplyr::select(features, target)

  train_rows = sample(nrow(ml_df)*train_ratio)
  train_df = ml_df[train_rows,]
  test_df = ml_df[-train_rows,]

  learner = makeLearner(model_type)
  pre_learner = makePreprocWrapperCaret(learner, method=c('center','scale'))

  train_task = makeRegrTask(data = train_df, target = target)
  test_task = makeRegrTask(data = test_df, target = target)
  all_task = makeRegrTask(data = ml_df, target = target)

  meas = list(mlr::mae, mlr::rmse)

  tune_res = tuneParams(
    learner, train_task, resampling=cv5, par.set=pset,
    control=makeTuneControlMBO(budget = mbo_budget),
    measures=meas
  )

  dir.create(paste0('../output/',prefix), showWarnings = FALSE)
  write_rds(tune_res, str_c('../output/',prefix,"/tune_results.rds"))

  # set hyperparameters
  tuned_learner = setHyperPars(learner = learner, par.vals = tune_res$x)

  # train final model for performance and interpretation
  model = mlr::train(tuned_learner, train_task)
  test_predict = predict(model, newdata = test_df)
  train_predict = predict(model, newdata = train_df)
  predictor = Predictor$new(model, data = ml_df)

  # monte-carlo performance measures
  resample_obj = mlr::resample(
    tuned_learner, all_task,
    makeResampleDesc("Subsample", iters=resample_iters, split=4/5, predict='both'),
    measures = meas,
    show.info = FALSE
  )

  model_res = get_resample_regr_res(resample_obj)
  model_res$train_perf = performance(predict(model, newdata=train_df), measures = meas)
  model_res$test_perf = performance(predict(model, newdata=test_df), measures = meas)
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

#' Make a 3D plot using Plot3D with a viridis colour scale and sf dataframe
#' @return plot
#' @export
make_3d_plot <- function(sf_df, z_col, col_col, col_func=viridis(50), size=0.1){
  coords = sf_df %>% sf::st_coordinates(.)
  scatter3D(
    x=coords[,1],
    y=coords[,2],
    z=sf_df %>% pull(z_col),
    colvar=sf_df %>% pull(col_col),
    pch=19, cex=size, bty='b2',
    theta=0, phi=30,
    ylab='N', xlab='E', zlab='',
    ticktype = 'detailed',
    col = col_func
  )
}

#' Wrapper of the Wrapper Function for model generation
#'
#' @param mldata the data (dataframe)
#' @param feats features for model (character vector)
#' @param feat_prefix the prefix for the feature set
#' @param targets named prefix of output directory
#' @return TRUE
#' @export
generate_models <- function(
  mldata,
  feats,
  feat_prefix,
  targets=c("youngs_modulus", "poisson_ratio", "brittleness")
){
  for (target in targets){
    if (target == 'youngs_modulus'){
      target_suffix='ym'
    } else if (target == 'poisson_ratio'){
      target_suffix='pr'
    } else {
      target_suffix='br'
    }

    # glm
    generate_model(
      data = mldata, target = target,
      features = feats, model_type = "regr.glmnet",
      pset = makeParamSet(
         makeNumericParam('alpha',lower = 0, upper = 1),
         makeIntegerParam('lambda', lower = -4, upper = 1, trafo = function(x) 10^x)
       ),
       prefix = str_c(feat_prefix,'_glm_',target_suffix)
    )

    # mars
    generate_model(
       data = mldata, target = target,
       features = feats, model_type = "regr.mars",
       pset = makeParamSet(
         makeDiscreteParam('degree',values=c(1,2,3)),
         makeIntegerParam('nk', lower = 1, upper = 10)
       ),
       prefix = str_c(feat_prefix,'_mars_',target_suffix)
    )

    # random forest
    generate_model(
       data = mldata, target = target,
       features = feats, model_type = "regr.ranger",
       pset= makeParamSet(
         makeIntegerParam('mtry', lower = 1L, upper = min(length(feats),6L)),
         makeIntegerParam('num.trees', lower = 10L, upper = 200L),
         makeIntegerParam('min.node.size', lower = 0, upper = 1, trafo = function(x) 80^x),
         makeNumericParam('sample.fraction', lower = 0.2, upper = 0.9)
       ),
       prefix = str_c(feat_prefix,'_rf_',target_suffix)
    )
  }
}

#' Wave Equation for Poisson Ratio with a correction factor
#'
#' @param df the data (dataframe)
#' @param correction_factor scalar correction factor (numeric)
#' @return vector, predicted poisson ratio values
#' @export
wave_pr_w_correction <- function(df, correction_factor){
  (df$pred_vp^2-2*df$pred_vs^2)/(2*(df$pred_vp^2-df$pred_vs^2)) - correction_factor
}

#' Wave Equation for Poisson Ratio with a correction factor
#'
#' @param df the data (dataframe)
#' @param correction_factor scalar correction factor (numeric)
#' @return vector, predicted poisson ratio values
#' @export
wave_ym_w_correction <- function(df, correction_factor){
  df$bulk_density*df$pred_vs^2*((3*df$pred_vp^2 - 4*df$pred_vs^2)/(df$pred_vp^2-df$pred_vs^2))/1E6 - correction_factor
}


#' Average of Mean Absolute Error and RMSE
#' of Wave Equation PR prediction with a correction factor
#' (to match statistical models)
#' Used for cross validated optimization
#'
#' @param df the data (dataframe)
#' @param correction_factor scalar correction factor (numeric)
#' @return mean absolute error
#' @export
pr_correction_loss <- function(correction_factor=0.01, df){
  rock_physics_pr = wave_pr_w_correction(df, correction_factor)
  mae = mean(abs(rock_physics_pr - df$poisson_ratio))
  rmse = sqrt(mean((rock_physics_pr - df$poisson_ratio)^2))
  mean(mae, rmse)
}

#' Average of Mean Absolute Error and RMSE
#' of Wave Equation YM prediction with a correction factor
#' (to match statistical models)
#' Used for cross validated optimization
#'
#' @param df the data (dataframe)
#' @param correction_factor scalar correction factor (numeric)
#' @return mean absolute error
#' @export
ym_correction_loss <- function(correction_factor=1, df){
  rock_physics_ym = wave_ym_w_correction(df, correction_factor)
  mae = mean(abs(rock_physics_ym - df$youngs_modulus))
  rmse = sqrt(mean((rock_physics_ym - df$youngs_modulus)^2))
  mean(mae, rmse)
}
