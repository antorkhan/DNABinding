featurefiltering <-
  function(features, rankedFeatures, maxFeatureCount = Inf) {
    columns = colnames(features)
    featureFilter = rankedFeatures[1:maxFeatureCount];
    featureFilter[length(featureFilter) + 1] = "protection"
    
    if (ncol(features) - 1 > maxFeatureCount) {
      for(i in 1:length(columns)){
        if(!columns[i] %in% featureFilter){
          features[columns[i]] = NULL
        }
      }
    }
    
    # Ensure that all the selected features do exist
    # in the input set (even if with 0 value)
    for(i in 1:length(featureFilter)){
      if(!featureFilter[i] %in% columns){
        features[featureFilter[i]] = 0
      }
    }
    
    return(features);
  }