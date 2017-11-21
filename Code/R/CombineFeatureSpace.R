library(e1071)
library(ROCR)
library(randomForest)

source('./combineFeatures.R');

timestamp();

fScheme = "_comb";

RDSFolder = "./"

rankedFeaturesFile = paste(RDSFolder, "ff"             , fScheme, ".rds", sep = "");
featureFile        = paste(RDSFolder, "featurized"     , fScheme, ".rds", sep = "");
testFeatureFile    = paste(RDSFolder, "testFeaturized" , fScheme, ".rds", sep = "");

fSubSchemes = c("_nGrams","_nGDip", "_PSF");

# Generate combined feature ranking
# if (!file.exists(rankedFeaturesFile)) {
#   combFeatureNames = c();
#   combFeatureScores = c();
#   for (fSubScheme in fSubSchemes) {
#     curRFmodelFile = paste("rfmodel_", as.character(nTrainData), fSubScheme, ".rds", sep = "");
#     curRFModel = readRDS(curRFmodelFile);
# 
#     curScores = scale(curRFModel$importance[,3]);
#     combFeatureScores = c(combFeatureScores, curScores);
#     combFeatureNames = c(combFeatureNames, rownames(curScores));
#   }
# 
#   saveRDS(combFeatureNames[order(-combFeatureScores)], rankedFeaturesFile);
# 
# } else {
#   cat(as.character(Sys.time()),">> Computing feature ranking ...\n");
#   rankedFeatures = readRDS(rankedFeaturesFile);
#   cat(as.character(Sys.time()),">> Done ( from cached file:", rankedFeaturesFile, ")\n");
# }

# Generate combined feature vectors for training and test sets

cat(as.character(Sys.time()),">> Combining Training set features ...\n");
if (!file.exists(featureFile)) {
  fileNamePrefix = paste0(RDSFolder,"featurized");
  features = combineFeatures(fSubSchemes, fileNamePrefix);
  saveRDS(features, featureFile);
  cat(as.character(Sys.time()),">> Done. Saved to ", featureFile, "\n");
} else {
  features = readRDS(featureFile);
  cat(as.character(Sys.time()),">> Done. (From cach file ", featureFile, ")\n");
}

cat(as.character(Sys.time()),">> Combining Test set features ...\n");
if (!file.exists(testFeatureFile)) {
  fileNamePrefix = paste0(RDSFolder, "testFeaturized");
  testFeatures = combineFeatures(fSubSchemes, fileNamePrefix);
  saveRDS(testFeatures, testFeatureFile);
  cat(as.character(Sys.time()),">> Done. Saved to ", testFeatureFile, "\n");
} else {
  testFeatures = readRDS(testFeatureFile);
  cat(as.character(Sys.time()),">> Done. (From cach file ", testFeatureFile, ")\n");
}