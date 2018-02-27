library(e1071)
library(ROCR)
library(pracma)
library(ggplot2)

source('./featurefiltering.R');
source('./learnWithCV.R');


timestamp();
set.seed(10);

# featureCountList = seq(from = 500, to = 6500, by = 500);
featureCountList = c(500); 

nFolds = 10

fScheme   = "_comb";

RDSFolder = "RDSFiles/"

# File names #
rankedFeaturesFile = paste(RDSFolder, "ff_SvmRFE" , fScheme, ".rds", sep = "");
featureFile        = paste(RDSFolder, "featurized", fScheme, ".rds", sep = "");
outFile            = paste("out"                  , fScheme, ".csv", sep = "");


cat(as.character(Sys.time()),">> Reading training set features from", featureFile, "...\n");
features = readRDS(featureFile);
cat(as.character(Sys.time()),">> Done\n");

cat(as.character(Sys.time()),">> Reading feature ranking from", rankedFeaturesFile, "...\n");
rankedFeatures = readRDS(rankedFeaturesFile);
cat(as.character(Sys.time()),">> Done\n");

# random shuffle of features
features <- features[sample(nrow(features)),]

bestPerf = NULL;
accData  = NULL;

cat(as.character(Sys.time()),">> Entering cross validation. Folds = ", nFolds, " ...\n");

# Reduce the feature vectors to the max size that we will be testing.
# This way the filtering cost in the loop below will be reduced.
features = featurefiltering(features, rankedFeatures, max(featureCountList));

# For regression study, we need to 'unfactor' the dependent var.
# When converting from factor to numeric, DNA-BPs becomes 2 and non DNA-BPs becomes 1.
# So we need to deduct 1.
features$protection = as.numeric(features$protection) - 1;

rocCurvePoints = NULL;
prCurvePoints = NULL;

for (maxFeatureCount in featureCountList) 
{
  trainingSet = featurefiltering(features, rankedFeatures, maxFeatureCount);

  perf = learnWithCV(protection ~ ., trainingSet, cross = nFolds, "svm", 
                     kernel = "linear", svmCost = svmC);
  
  AUCROC = perf$AUCROC;
  df = data.frame(
        x = unlist(perf$rocCurve@x.values), 
        y = unlist(perf$rocCurve@y.values), 
        Features = as.character(maxFeatureCount)
        );
  rocCurvePoints = rbind(rocCurvePoints, df);
  
  cat(maxFeatureCount, ",", AUCROC, "\n");
  accData = rbind(accData, c(maxFeatureCount, AUCROC));
  colnames(accData) = c("nF", "AUCROC");

  write.csv(accData, outFile);
    
  cat("\n");
}

rocCurveFile = "ROCCurve.eps";

rocPlot = ggplot(rocCurvePoints,aes(x, y)) + 
  theme_bw(base_size = 36, base_family = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top") +
  #theme(aspect.ratio = 0.7) +
  geom_line(aes(colour=Features),size = 3) +
  labs(x = "False Positive Rate", y = "True Positive Rate");
  
postscript(file = rocCurveFile, paper = "letter");
rocPlot;
dev.off();
