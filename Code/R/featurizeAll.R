library(e1071)
library(ROCR)
library(randomForest)

source('./featurization.R');

timestamp();

dataFile          = "trainingSet.csv";
featureFilePrefix = "featurized";
fScheme           = "_nGrams";

featureFile = paste(featureFilePrefix, fScheme, ".rds", sep = "");

amins = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y");

cat(as.character(Sys.time()),">> Featurizing ...\n");
if (!file.exists(featureFile)) {
  cat(as.character(Sys.time()),">> Reading sequence file (", dataFile, ") ...\n");
  data = readRDS(dataFile);
  nData = length(data[,1]);
  cat(as.character(Sys.time()),">> Done\n");
  
  features = featurization(data$Sequence, data$protection, amins, nGramOrder = 1, nGDipOrder = 0, psfOrder = 0);
  
  features = cbind("ID" = data$ID, features);
  saveRDS(features, featureFile);
  cat(as.character(Sys.time()),">> Featurizing Done.\n");
} else {
  features = readRDS(featureFile);
  cat(as.character(Sys.time()),">> Done ( from cached file:", featureFile, ")\n");
}
cat(as.character(Sys.time()),">> Total features: ", length(features[1,]) - 2, "\n");
