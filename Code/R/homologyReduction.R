homologyReduction <-
  function(features, hrScheme) {
    result = NULL;
    if (is.null(hrScheme) || hrScheme == "") {
      result = features;
    } else {
      homologyRemovedFile = paste("trainingSet" , hrScheme, ".csv", sep = "");
      novelSet = read.csv(homologyRemovedFile);
    
      result = features[which (features$ID %in% novelSet$ID),];
    }
    
    return(result);
  }