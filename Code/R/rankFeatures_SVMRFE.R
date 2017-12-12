library(e1071)

source("svmRFE.R")
source("featurefiltering.R")
source("homologyReduction.R")

timestamp();

set.seed(10);

fScheme  = "_comb";
hrScheme = "_BLASTCLUST25"
bScheme  = "_RUS";

RDSFolder = "./"

fileNameSuffix = paste(fScheme, ".rds", sep = "");

InitialRankedFeaturesFile = paste(RDSFolder, "ff"        , hrScheme,          fScheme, ".rds", sep = "");
FinalRankedFeaturesFile   = paste(RDSFolder, "ff_SvmRFE" , hrScheme, bScheme, fScheme, ".rds", sep = "");
featureFile               = paste(RDSFolder, "featurized",                    fScheme, ".rds", sep = "");

if (!file.exists(FinalRankedFeaturesFile)) {
  cat(as.character(Sys.time()),">> Loading feature file ...\n");
  features = readRDS(featureFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", featureFile, ")\n");
  
  cat(as.character(Sys.time()),">> Removing homology. hrScheme = ", hrScheme, "...\n");
  features = homologyReduction(features, hrScheme);
  cat(as.character(Sys.time()),">> Done\n");
  
  features$ID = NULL;
  cat(as.character(Sys.time()),">> Total features: ", length(features[1,]) - 1, "\n");

  cat(as.character(Sys.time()),">> Loading initial feature ranking ...\n");
  rankedFeatures = readRDS(InitialRankedFeaturesFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", InitialRankedFeaturesFile, ")\n");
  
  # Reduce the feature vectors to the max size that we will be testing.
  features = featurefiltering(features, rankedFeatures, 7000);
  
  cat(as.character(Sys.time()),">> Balancing scheme is:", bScheme, "\n");
  if (bScheme != "") {
    cat(as.character(Sys.time()),">> Balancing ...\n");
    
    #
   # Balance the training set by undersampling the larger set
    #
    nPositive = length(which(features$protection == 1));
    nNegative = length(features[,1]) - nPositive;
    nBalanced = min(nPositive, nNegative);
    
    positiveSetInd = sample(1:nPositive)[1:nBalanced];
    negativeSetInd = sample((nPositive+1):length(features[,1]))[1:nBalanced];
    
    features = rbind(features[positiveSetInd,], features[negativeSetInd,]);

    cat(as.character(Sys.time()),">> Done.\n");
  }
  
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
