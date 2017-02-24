##################################################################
##################################################################
#####                                                         ####
#####              Super Learner parameters                   ####
#####                                                         ####
##################################################################
##################################################################


### xgboost component for Super Learner

set.seed(1)
SL.xgboost = function(Y, X, newX, family, obsWeights, id, ntrees = 180, #100,
                      max_depth = 4,
                      shrinkage = 0.1, minobspernode = 10,
                      params = list(),
                      nthread = 1, verbose = 0,
                      ...) {
    
  # Convert to an xgboost compatible data matrix, using the sample weights.
  xgmat = xgboost::xgb.DMatrix(data=as.matrix(X), label=Y, weight = obsWeights)

  # TODO: support early stopping, which requires a "watchlist". See ?xgb.train

  if (family$family == "gaussian") {
    model = xgboost::xgboost(data=xgmat, objective="reg:linear", nround = ntrees,
                             max_depth = max_depth, minchildweight = minobspernode, eta = shrinkage, verbose=verbose,
                             nthread = nthread, params = params)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data=xgmat, objective="binary:logistic", nround = ntrees,
                             max_depth = max_depth, minchildweight = minobspernode, eta = shrinkage, verbose=verbose,
                             nthread = nthread, params = params)
  }
  if (family$family == "multinomial") {
    # TODO: test this.
    model = xgboost::xgboost(data=xgmat, objective="multi:softmax", nround = ntrees,
                             max_depth = max_depth, minchildweight = minobspernode, eta = shrinkage, verbose=verbose,
                             num_class=length(unique(Y)), nthread = nthread, params = params)
  }
  pred = predict(model, newdata=data.matrix(newX))
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

predict.SL.xgboost = function (object, newdata, ...) {
  pred = predict(object$object, newdata = data.matrix(newdata))
  return(pred)
}
