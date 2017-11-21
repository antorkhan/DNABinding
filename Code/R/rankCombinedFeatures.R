library(randomForest)

source("featurefiltering.R")

timestamp();

fScheme = "_comb";

RDSFolder = "./"

rfmodelFile        = paste(RDSFolder, "rfmodel"        , fScheme, ".rds", sep = "");
rankedFeaturesFile = paste(RDSFolder, "ff"             , fScheme, ".rds", sep = "");
featureFile        = paste(RDSFolder, "featurized"     , fScheme, ".rds", sep = "");

fSubSchemes = c("_nGrams", "_nGDip", "_PSF");

# Generate combined feature ranking
if (!file.exists(rankedFeaturesFile)) {
  
  if (!file.exists(rfmodelFile)) {
    combFeatureNames = c();
    
    cat(as.character(Sys.time()),">> Heuristically combining feature space(s) ...\n");
    for (fSubScheme in fSubSchemes) {
      curRFmodelFile = paste("rfmodel", fSubScheme, ".rds", sep = "");
      curRFModel = readRDS(curRFmodelFile);
      
      sortedImpScore = curRFModel$importance[order(-curRFModel$importance[,3]), 3];
      lastSelected = which(sortedImpScore <= 0)[1] - 1;
      
      names(sortedImpScore[1:lastSelected]);
      combFeatureNames = c(combFeatureNames, names(sortedImpScore[1:lastSelected]));
    }
    cat(as.character(Sys.time()),">> Done.\n");
    
    cat(as.character(Sys.time()),">> Loading feature file ...\n");
    features = readRDS(featureFile);
    cat(as.character(Sys.time()),">> Done. ( from cached file:", featureFile, ")\n");
    
    # Reduce the feature vector to the heuristically combined feature space
    cat(as.character(Sys.time()),">> Keeping heuristically combined features ...\n");
    features = featurefiltering(features, combFeatureNames, length(combFeatureNames));
    cat(as.character(Sys.time()),">> Done. Total features: ", length(features[1,]) - 1, "\n");
    
    cat(as.character(Sys.time()),">> Computing random forest ...\n");

    rfmodel = randomForest(protection ~ ., features, importance=TRUE);
    saveRDS(rfmodel, rfmodelFile);
    cat(as.character(Sys.time()),">> Done.\n");
  } else {
    rfmodel = readRDS(rfmodelFile);
    cat(as.character(Sys.time()),">> Done ( from cached file:", rfmodelFile, ")\n");
  }
  
  cat(as.character(Sys.time()),">> Computing feature ranking ...\n");
  rankedFeatures = rownames(rfmodel$importance[order(-rfmodel$importance[,3]),])
  saveRDS(rankedFeatures, rankedFeaturesFile);
  cat(as.character(Sys.time()),">> Done\n");
} else {
  cat(as.character(Sys.time()),">> Computing feature ranking ...\n");
  rankedFeatures = readRDS(rankedFeaturesFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", rankedFeaturesFile, ")\n");
}