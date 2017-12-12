library(randomForest)

source("./homologyReduction.R")

timestamp();

set.seed(10);

fScheme = "_nGrams";
hrScheme = "_BLASTCLUST25"

RDSFolder          = "RDSFiles/"

rfmodelFile        = paste(RDSFolder, "rfmodel"   , hrScheme, fScheme, ".rds", sep = "");
rankedFeaturesFile = paste(RDSFolder, "ff"        , hrScheme, fScheme, ".rds", sep = "");
featureFile        = paste(RDSFolder, "featurized",           fScheme, ".rds", sep = "");


if (!file.exists(rankedFeaturesFile)) {
  cat(as.character(Sys.time()),">> Loading feature file ...\n");
  features = readRDS(featureFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", featureFile, ")\n");
  
  cat(as.character(Sys.time()),">> Removing homology. hrScheme = ", hrScheme, "...\n");
  features = homologyReduction(features, hrScheme);
  cat(as.character(Sys.time()),">> Done\n");
  
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
