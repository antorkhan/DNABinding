library(e1071)

source("svmRFE.R")
source("featurefiltering.R")

timestamp();

set.seed(10);

fScheme = "_comb";

RDSFolder = "RDSFiles/"

fileNameSuffix = paste(fScheme, ".rds", sep = "");

InitialRankedFeaturesFile = paste(RDSFolder, "ff" , fileNameSuffix, sep = "");
FinalRankedFeaturesFile   = paste(RDSFolder, "ff_SvmRFE", fileNameSuffix, sep = "");
featureFile               = paste(RDSFolder, "featurized", fileNameSuffix, sep = "");

if (!file.exists(FinalRankedFeaturesFile)) {
  cat(as.character(Sys.time()),">> Loading feature file ...\n");
  features = readRDS(featureFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", featureFile, ")\n");
 
  features$ID = NULL;
  cat(as.character(Sys.time()),">> Total features: ", length(features[1,]) - 1, "\n");

  cat(as.character(Sys.time()),">> Loading initial feature ranking ...\n");
  rankedFeatures = readRDS(InitialRankedFeaturesFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", InitialRankedFeaturesFile, ")\n");
  
  # Reduce the feature vectors to the max size that we will be testing.
  # This way the filtering cost in the loop below will be reduced.
  features = featurefiltering(features, rankedFeatures, 7000);
  
  # random shuffle of features
  features <- features[sample(nrow(features)),]
  
  cat(as.character(Sys.time()),">> Computing feature ranking ...\n");
  
  labelCol = which(colnames(features) == "protection");
  
  rankedFeatures = svmRFE(features[,-labelCol], features$protection, 25);
  saveRDS(rankedFeatures, FinalRankedFeaturesFile);
  cat(as.character(Sys.time()),">> Done\n");
} else {
  cat(as.character(Sys.time()),">> Computing feature ranking ...\n");
  rankedFeatures = readRDS(FinalRankedFeaturesFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", FinalRankedFeaturesFile, ")\n");
}