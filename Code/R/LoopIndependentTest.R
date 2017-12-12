library(e1071)
library(ROCR)

source('featurefiltering.R');
source('learn.R');

timestamp();

set.seed(10);

fScheme = "_comb";
hrScheme = "_BLASTCLUST25";

featureCountList = seq(from=100, to=1500, by=100);
#featureCountList = c(200);

# File names #
outFile     = "IndependentTestResults.csv";

RDSFolder          = "RDSFiles/"
rankedFeaturesFile = paste(RDSFolder, "ff_SvmRFE"    , hrScheme, fScheme, ".rds", sep = "");
featureFile        = paste(RDSFolder, "featurized"    ,           fScheme, ".rds", sep = "");
testFeatureFile    = paste(RDSFolder, "testFeaturized",           fScheme, ".rds", sep = "");

cat(as.character(Sys.time()),">> Reading feature ranking from", rankedFeaturesFile, "...\n");
rankedFeatures = readRDS(rankedFeaturesFile);
cat(as.character(Sys.time()),">> Done\n");

cat(as.character(Sys.time()),">> Reading training set features from", featureFile, "...\n");
features = readRDS(featureFile);
cat(as.character(Sys.time()),">> Done\n");

cat(as.character(Sys.time()),">> Removing homology. hrScheme = ", hrScheme, "...\n");
features = homologyReduction(features, hrScheme);
cat(as.character(Sys.time()),">> Done\n");

# random shuffle of features
features <- features[sample(nrow(features)),]

cat(as.character(Sys.time()),">> Reading test set features from", testFeatureFile, "...\n");
testSet = readRDS(testFeatureFile);
cat(as.character(Sys.time()),">> Done\n");

bestPerf = NULL;
bestParams = NULL;
accData = NULL;

cat(as.character(Sys.time()),">> Entering independet testing ...\n");

# Reduce the feature vectors to the max size that we will be testing.
# This way the filtering cost in the loop below will be reduced.
features = featurefiltering(features, rankedFeatures, max(featureCountList));
testSet  = featurefiltering(testSet, rankedFeatures, max(featureCountList));

# For regression study, we need to 'unfactor' the dependent var.
# When converting from factor to numeric, Antigens becomes 2 and Non-antigens becomes 1.
# So we need to deduct 1.
features$protection = as.numeric(features$protection) - 1;
testSet$protection  = as.numeric(testSet$protection) - 1;

for (maxFeatureCount in featureCountList) 
{
  trainingSet = featurefiltering(features, rankedFeatures, maxFeatureCount);
  
  rfModel = learn(protection ~ ., trainingSet, "svm", kernel = "linear");
  rfPred = predict(rfModel, testSet);
  rfPrediction = prediction(rfPred, testSet$protection);
  
  # Find the AUCROC
  auc  = ROCR::performance(rfPrediction,"auc")@y.values[[1]];
  
  threshold = 0.5;
  rfPrediction = prediction(as.numeric(rfPred >= threshold), testSet$protection);
  acc  = unlist(ROCR::performance(rfPrediction,"acc")@y.values)[2];
  sens = unlist(ROCR::performance(rfPrediction,"sens")@y.values)[2];
  spec = unlist(ROCR::performance(rfPrediction,"spec")@y.values)[2];
  mcc  = unlist(ROCR::performance(rfPrediction,"mat")@y.values)[2];
  prec = unlist(ROCR::performance(rfPrediction,"prec")@y.values)[2];
  f1   = unlist(ROCR::performance(rfPrediction,"f")@y.values)[2];
  
  cat(
    maxFeatureCount, 
    ",", round(auc, 2),
    ",", round(acc, 2),
    ",", round(sens, 2),
    ",", round(spec, 2),
    ",", round(prec, 2),
    ",", round(f1, 2),
    ",", round(mcc, 2),
    "\n"
  );
  
  accData = rbind(accData, c(maxFeatureCount, auc, acc, sens, spec, prec, f1, mcc));
  colnames(accData) = c("nF", "AUCROC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1" , "MCC");
  write.csv(accData, outFile);
}