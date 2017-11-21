library(randomForest)

timestamp();

set.seed(10);

fScheme = "_nGrams";

fileNameSuffix = paste(fScheme, ".rds", sep = "");

rfmodelFile        = paste("rfmodel"   , fileNameSuffix, sep = "");
rankedFeaturesFile = paste("ff"        , fileNameSuffix, sep = "");
featureFile        = paste("featurized", fileNameSuffix, sep = "");


if (!file.exists(rankedFeaturesFile)) {
  cat(as.character(Sys.time()),">> Loading feature file ...\n");
  features = readRDS(featureFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", featureFile, ")\n");
 
  features$ID = NULL;
  cat(as.character(Sys.time()),">> Total features: ", length(features[1,]) - 1, "\n");
  
  cat(as.character(Sys.time()),">> Computing random forest ...\n");
  if (!file.exists(rfmodelFile)) {
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
