require("e1071");
require("randomForest");
require("ROCR");
require("pracma");

source("learn.R");

learnWithCV <-
  function(formula, data, cross, learner, ...) {
    N = length(data[, 1])
    folds = seq(from=1,to=N, by=round(N/cross))
    folds[cross+1] = N+1
    predVector = c()
    for (i in 1:cross) {
      trainFolds = data
      testFold = data[(folds[i]:(folds[i+1]-1)),]
      trainFolds = data[-(folds[i]:(folds[i+1]-1)),]

      model = learn(formula, trainFolds, learner, ...);
      mlPred = predict(model, testFold)
      predVector = c(predVector, as.numeric(mlPred))
      i = i + 1
    }
    
    # Now generate the model on full dataset to find the no. of support vectors
    # This can be used for model selection in. Lesser number of SVs will
    # result in better generalization
    model = learn(formula, data, learner, ...);

    dependentVar = all.vars(formula)[1];
    
    # perform classification based perf. measures
    if (is.factor(data[,dependentVar])) {
      mlPrediction = prediction(as.numeric(predVector), as.numeric(data[,dependentVar]))
      
      acc = unlist(ROCR::performance(mlPrediction,"acc")@y.values)[2]
      sensitivity = unlist(ROCR::performance(mlPrediction,"sens")@y.values)[2];
      specificity = unlist(ROCR::performance(mlPrediction,"spec")@y.values)[2];
      precision   = unlist(ROCR::performance(mlPrediction,"prec")@y.values)[2];
      mcc = unlist(ROCR::performance(mlPrediction,"mat")@y.values)[2];
      f1  = unlist(ROCR::performance(mlPrediction,"f")@y.values)[2];
      
      return(list(
        "model"= model,
        "acc"  = acc,
        "sens" = sensitivity,
        "spec" = specificity,
        "mcc"  = mcc,
        "prec" = precision,
        "f1"   = f1
      ))
    }
    else {
      # perform regression based perf. measurements
      mlPrediction = prediction(predVector, data[, dependentVar]);
      
      # Find the ROC curve and AUCROC
      AUCROC  = ROCR::performance(mlPrediction,"auc")@y.values[[1]];
      rocCurve = ROCR::performance(mlPrediction,"tpr", "fpr");

      # Find the PR curve and AUCPR
      prCurve  = ROCR::performance(mlPrediction,"prec", "rec");
      x = unlist(prCurve@x.values);
      y = unlist(prCurve@y.values);
      df = data.frame(x = x[2:length(x)], y = y[2:length(y)]);
      AUCPR  = trapz(df$x, df$y)
      
      # Find optimal threshold based on accuracy
      # accSeries = ROCR::performance(mlPrediction,"acc");
      # threshold = unlist(accSeries@x.values)[[which.max(unlist(accSeries@y.values))]];
      
      # Use a fixed threshold of 0.5
      threshold = 0.5;
      
      mlPrediction = prediction(as.numeric(predVector >= threshold), data[, dependentVar]);
      
      acc = unlist(ROCR::performance(mlPrediction,"acc")@y.values)[2]
      sensitivity = unlist(ROCR::performance(mlPrediction,"sens")@y.values)[2];
      specificity = unlist(ROCR::performance(mlPrediction,"spec")@y.values)[2];
      precision   = unlist(ROCR::performance(mlPrediction,"prec")@y.values)[2];
      mcc = unlist(ROCR::performance(mlPrediction,"mat")@y.values)[2];
      f1  = unlist(ROCR::performance(mlPrediction,"f")@y.values)[2];
      
      return(list(
        "model"     = model,
        "threshold" = threshold,
        "rocCurve"  = rocCurve,
        "prCurve"   = prCurve,
        "AUCROC"    = AUCROC,
        "AUCPR"     = AUCPR,
        "acc"       = acc,
        "sens"      = sensitivity,
        "spec"      = specificity,
        "prec"      = precision,
        "f1"        = f1,
        "mcc"       = mcc
      ))
    }
  }