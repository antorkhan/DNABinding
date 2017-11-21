combineFeatures <-
  function(featureSchemes, subFeatureFilePrefix) {
    features = NULL;
    for (fScheme in featureSchemes) {
      curFeatureFile = paste(subFeatureFilePrefix, fScheme, ".rds", sep = "");
      curFeatures = readRDS(curFeatureFile);
      curFeatures$Serial <- seq.int(nrow(curFeatures))
      cat(as.character(Sys.time()),">> Read file: ", curFeatureFile, "\n");
      # ID and protection are common columns in each sub feature file
      # We will merge by ID. The remaining 1 column should come from the
      # first file. Subsequently they will be ignored from other files
      if (is.null(features)) {
        features = curFeatures;
      } else {
        curFeatures$protection = NULL;
        curFeatures$ID = NULL;
        features = merge(features, curFeatures, by="Serial");
      }
    }
    features$Serial = NULL;
    cat(as.character(Sys.time()),">> Total features: ", length(features[1,]) - 2, "\n");
    return(features)
  }